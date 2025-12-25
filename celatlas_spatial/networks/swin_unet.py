import os
import torch

from celatlas_spatial.networks.vision_transformer import SwinUnet as ViT_seg

class SwinUnetImpl(object):
    def __init__(self, args, model_path=''):
        self._device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self._model_path = model_path
        self._model = ViT_seg(args.config, img_size=args.img_size, num_classes=args.num_classes).to(self._device)
        self._f_load_model()

    def _f_load_model(self):
        if os.path.exists(self._model_path):
            state_dict = torch.load(self._model_path, map_location=self._device)
            self._model.load_state_dict(state_dict, strict=False)

    def f_predict(self, img):
        self._model.eval()
        with torch.no_grad():
            output = self._model(img)
            output = torch.argmax(torch.softmax(output, dim=1), dim=1).squeeze(0)
            output = output.cpu().detach().numpy()
        return output
