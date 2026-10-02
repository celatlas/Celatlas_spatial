import cv2
import torch
import numpy as np

from scipy.ndimage import zoom
from celatlas_spatial.networks import swin_unet
from celatlas_spatial.tools.config import get_config


class Args:
    def __init__(self, img_size, num_classes):
        self.img_size = img_size
        self.num_classes = num_classes
        self.config = get_config()

class SwinChipCut(object):
    def __init__(self, image_size=224, model_path=''):
        self._device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self._args = Args(image_size, 2)
        self._image_size = image_size
        self._model_path = model_path
        self._sw_model = None

        self._init_model()

    def _init_model(self):
        self._sw_model = swin_unet.SwinUnetImpl(self._args, self._model_path)

    def f_predict(self, img):
        if len(img.shape) == 2:
            x, y = img.shape
            if x != self._image_size or y != self._image_size:
                img = zoom(img, (self._image_size / x, self._image_size / y), order=3)
            img = cv2.cvtColor(img, cv2.COLOR_GRAY2RGB).transpose((2, 0, 1))
            img = torch.from_numpy(img).unsqueeze(0).float().to(self._device)
            pred = self._sw_model.f_predict(img)
            pred = zoom(pred, (x / self._image_size, y / self._image_size), order=3)
        elif len(img.shape) == 3:
            x, y, _ = img.shape
            if x != self._image_size or y != self._image_size:
                img = zoom(img, (self._image_size / x, self._image_size / y, 1), order=3)
            img = img.transpose((2, 0, 1))
            img = torch.from_numpy(img).unsqueeze(0).float().to(self._device)
            pred = self._sw_model.f_predict(img)
            pred = zoom(pred, (x / self._image_size, y / self._image_size), order=3)
        else:
            raise ValueError("Error: image shape is not as expected！")
        return pred
