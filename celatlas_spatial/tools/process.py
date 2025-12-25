import cv2
import tifffile
import numpy as np

# Set matplotlib to non-interactive backend before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class Processor:
    def __init__(self, output, image):
        self.output = output
        self.image = image

    def rotate_correction(self):
        pass

    def detect_border(self, border_value=179, tolerance=5, min_area=100, debug=False):
        gray = self.image.copy() if len(self.image.shape) != 3 else cv2.cvtColor(self.image, cv2.COLOR_BGR2GRAY)

        mask = np.zeros(gray.shape, dtype=np.uint8)
        binary = cv2.inRange(gray, border_value - tolerance, border_value + tolerance)
        contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        detected_shapes = []
        contour = max(contours, key=cv2.contourArea)
        area = cv2.contourArea(contour)
        if area > min_area:
            perimeter = cv2.arcLength(contour, True)

            epsilon = 0.02 * perimeter
            approx = cv2.approxPolyDP(contour, epsilon, True)

            if len(approx) == 4:
                cv2.fillPoly(mask, [approx], 255)
                detected_shapes = [approx]
            else:
                print(f"Warning: Detected shape has {len(approx)} vertices, expected 4.")

        if debug:
            print(f"Image shape: {self.image.shape}")
            print(f"Number of detected shapes: {len(detected_shapes)}")

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

            # Original image with detected shapes
            ax1.imshow(cv2.cvtColor(self.image, cv2.COLOR_BGR2RGB))
            for shape in detected_shapes:
                ax1.add_patch(plt.Polygon(shape.reshape(-1, 2), fill=False, edgecolor='r', linewidth=2))
            ax1.set_title('Original Image with Detected Shapes')

            # Mask image
            ax2.imshow(mask, cmap='gray')
            ax2.set_title('Mask')

            plt.show()

            # Print information about each detected shape
            for i, shape in enumerate(detected_shapes):
                print(f"Shape {i + 1}:")
                print(f"  Number of vertices: {len(shape)}")
                print(f"  Area: {cv2.contourArea(shape)}")
                print(f"  Perimeter: {cv2.arcLength(shape, True)}")

                # Calculate and print bounding box information
                rect = cv2.minAreaRect(shape)
                box = cv2.boxPoints(rect)
                box = np.intp(box)
                print(f"  Bounding box:")
                print(f"    Center: {rect[0]}")
                print(f"    Size: {rect[1]}")
                print(f"    Angle: {rect[2]}")
                print()
        return mask

    def detect_polygons(self, border_value, tolerance=5, min_area=100):
        gray = self.image.copy() if len(self.image.shape) != 3 else cv2.cvtColor(self.image, cv2.COLOR_BGR2GRAY)
        mask = np.zeros(gray.shape, dtype=np.uint8)
        binary = cv2.inRange(gray, border_value - tolerance, 255)
        contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        detected_polygons = []

        contour = max(contours, key=cv2.contourArea)
        area = cv2.contourArea(contour)
        if area > min_area:
            perimeter = cv2.arcLength(contour, True)
            epsilon = 0.02 * perimeter
            approx = cv2.approxPolyDP(contour, epsilon, True)  # approximate the contour

            vertices = len(approx)
            if vertices >= 3:
                detected_polygons.append({
                    'contour': approx,
                    'area': area,
                })

            cv2.drawContours(mask, [approx], 0, 255, -1)
        return mask, detected_polygons

    @staticmethod
    def get_nonzero_bbox(mask):
        nonzero = cv2.findNonZero(mask)
        x, y, w, h = cv2.boundingRect(nonzero)
        return x, y, w, h

    @staticmethod
    def get_region_center(image):
        M = cv2.moments(image)
        cx = int(M['m10'] / M['m00'])
        cy = int(M['m01'] / M['m00'])
        return cx, cy

    @staticmethod
    def get_rot_point(point, angle, center):
        x, y = point[0] - center[0], point[1] - center[1]
        rho = np.sqrt(x ** 2 + y ** 2)
        theta = np.arctan2(y, x)

        new_theta = theta + np.radians(angle)
        new_x = rho * np.cos(new_theta) + center[0]
        new_y = rho * np.sin(new_theta) + center[1]

        return np.array([new_x, new_y]).astype(np.float32)

    @staticmethod
    def rotate_points(points, angle, center):
        points = points - center
        angle_rad = np.radians(angle)

        cos_theta, sin_theta = np.cos(angle_rad), np.sin(angle_rad)
        rotation_matrix = np.array([[cos_theta, -sin_theta],
                                    [sin_theta, cos_theta]])

        rotated_points = np.dot(points, rotation_matrix.T)
        return rotated_points + center

    @staticmethod
    def get_boundary_point(image):
        mask = image.copy()
        mask[image != 0] = 255
        mask = cv2.dilate(mask, np.ones((3, 3), np.uint8), iterations=5)

        boundary_points = []
        contours, hierarchy = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        contours = sorted(contours, key=cv2.contourArea, reverse=True)[:4]
        for contour in contours:
            # use Douglas-Peucker algorithm to approximate the contour
            epsilon = 0.009 * cv2.arcLength(contour, True)
            approx = cv2.approxPolyDP(contour, epsilon, True)
            # if the approximated contour has four(or more) points, then assume that screen is found
            # make sure it is a square at least
            if len(approx) >= 4:
                boundary_points.extend(approx)

        return boundary_points

    def get_offset(self, image):
        boundary_points = self.get_boundary_point(image)
        boundary_points = np.array(boundary_points)
        left_top_index = np.linalg.norm(boundary_points - np.array([0, 0]), axis=2).argmin()
        transform_dis = (boundary_points[left_top_index].squeeze() - 100).clip(min=0)
        return transform_dis.astype(np.int32)

    def find_top_line(self, image):
        # find the boundary of the tissue
        boundary_points = self.get_boundary_point(image)
        boundary_points.sort(key=lambda point: point[0][1])  # order by y coordinate
        center = self.get_region_center(image)
        point1, point2 = boundary_points[:2]  # get the first two points

        # calculate the slope of the line
        x1, y1 = point2[0]  # top left point
        x2, y2 = point1[0]  # top right point
        if x2 - x1 != 0:
            slope = (y2 - y1) / (x2 - x1)
            angle_rad = np.arctan(slope)
            angle_deg = np.degrees(angle_rad)
        else:
            slope = float('inf')
            angle_rad = np.pi / 2
            angle_deg = 90.0
        return angle_deg, slope, center

    def find_bottom_line(self, image):
        # find the boundary of the tissue
        boundary_points = self.get_boundary_point(image)
        boundary_points.sort(key=lambda point: point[0][1])  # order by y coordinate
        center = self.get_region_center(image)
        point1, point2 = boundary_points[-2:]  # get the last two points

        # calculate the slope of the line
        x1, y1 = point1[0]
        x2, y2 = point2[0]
        if x2 - x1 != 0:
            slope = (y2 - y1) / (x2 - x1)
            angle_rad = np.arctan(slope)
            angle_deg = np.degrees(angle_rad)
        else:
            slope = float('inf')
            angle_rad = np.pi / 2
            angle_deg = 90.0
        return angle_deg, slope, center

    def get_rotate_image(self, image, angle, center):
        h, w = image.shape[:2]
        cos = np.abs(np.cos(np.radians(angle)))
        sin = np.abs(np.sin(np.radians(angle)))
        new_w = int((h * sin) + (w * cos))
        new_h = int((h * cos) + (w * sin))

        tx = (new_w - w) / 2
        ty = (new_h - h) / 2

        new_center = (center[0] + tx, center[1] + ty)

        rotation = cv2.getRotationMatrix2D(new_center, angle, 1.0)
        rotated = cv2.warpAffine(image, rotation, (new_w, new_h))
        return rotated, new_center

    def write_rotate_image(self, angle, center):
        h, w = self.image.shape[:2]
        cos = np.abs(np.cos(np.radians(angle)))
        sin = np.abs(np.sin(np.radians(angle)))
        new_w = int((h * sin) + (w * cos))
        new_h = int((h * cos) + (w * sin))

        tx = (new_w - w) / 2
        ty = (new_h - h) / 2

        new_center = (center[0] + tx, center[1] + ty)

        rotation = cv2.getRotationMatrix2D(new_center, angle, 1.0)
        rotated = cv2.warpAffine(self.image, rotation, (new_w, new_h))
        tifffile.imwrite(self.output, rotated, photometric="minisblack")
