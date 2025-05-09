import torch
import torch.nn as nn
import torch.nn.functional as F

class Conv2dWithPadding(nn.Conv2d):
    def __init__(self, *args, padding_mode='constant', padding_value=1, **kwargs):
        super(Conv2dWithPadding, self).__init__(*args, **kwargs)
        self.padding_mode = padding_mode

    def forward(self, x):
        if self.padding_mode == 'constant':
            x = F.pad(x, (1, 1, 1, 1), mode=self.padding_mode, value=self.padding_value)
        else:
            x = F.pad(x, (1, 1, 1, 1), mode=self.padding_mode)
        return super(Conv2dWithPadding, self).forward(x)

class GeneratorUNet(nn.Module):
    def __init__(self, in_channels=1, out_channels=1, init_features=16):
        super(GeneratorUNet, self).__init__()

        features = init_features
        self.encoder1 = self._block(in_channels, features, name="enc1")
        self.pool1 = nn.MaxPool2d(kernel_size=2, stride=2)
        self.encoder2 = self._block(features, features * 2, name="enc2")
        self.pool2 = nn.MaxPool2d(kernel_size=2, stride=2)
        self.encoder3 = self._block(features * 2, features * 4, name="enc3")
        self.pool3 = nn.MaxPool2d(kernel_size=2, stride=2)
        self.encoder4 = self._block(features * 4, features * 8, name="enc4")
        self.pool4 = nn.MaxPool2d(kernel_size=2, stride=2)

        self.bottleneck = self._block(features * 8, features * 16, name="bottleneck")

        self.upconv4 = nn.ConvTranspose2d(features * 16, features * 8, kernel_size=2, stride=2)
        self.decoder4 = self._block((features * 8) * 2, features * 8, name="dec4")
        self.upconv3 = nn.ConvTranspose2d(features * 8, features * 4, kernel_size=2, stride=2)
        self.decoder3 = self._block((features * 4) * 2, features * 4, name="dec3")
        self.upconv2 = nn.ConvTranspose2d(features * 4, features * 2, kernel_size=2, stride=2)
        self.decoder2 = self._block((features * 2) * 2, features * 2, name="dec2")
        self.upconv1 = nn.ConvTranspose2d(features * 2, features, kernel_size=2, stride=2)
        self.decoder1 = self._block(features * 2, features, name="dec1")

        self.conv = nn.Conv2d(in_channels=features, out_channels=out_channels, kernel_size=1)

    def forward(self, x):
        # x = x.unsqueeze(dim=1)
        enc1 = self.encoder1(x)
        enc2 = self.encoder2(self.pool1(enc1))
        enc3 = self.encoder3(self.pool2(enc2))
        enc4 = self.encoder4(self.pool3(enc3))

        bottleneck = self.bottleneck(self.pool4(enc4))

        dec4 = self.upconv4(bottleneck)
        dec4 = torch.cat((dec4, enc4), dim=1)
        dec4 = self.decoder4(dec4)
        dec3 = self.upconv3(dec4)
        dec3 = torch.cat((dec3, enc3), dim=1)
        dec3 = self.decoder3(dec3)
        dec2 = self.upconv2(dec3)
        dec2 = torch.cat((dec2, enc2), dim=1)
        dec2 = self.decoder2(dec2)
        dec1 = self.upconv1(dec2)
        dec1 = torch.cat((dec1, enc1), dim=1)
        dec1 = self.decoder1(dec1)

        return torch.sigmoid(self.conv(dec1))

    @staticmethod
    def _block(in_channels, features, name):
        return nn.Sequential(
            Conv2dWithPadding(
                in_channels=in_channels,
                out_channels=features,
                kernel_size=3,
                bias=False,
                padding_mode='replicate'
            ),
            nn.BatchNorm2d(num_features=features),
            nn.ReLU(inplace=False),
            Conv2dWithPadding(
                in_channels=features,
                out_channels=features,
                kernel_size=3,
                bias=False,
                padding_mode='replicate'
            ),
            nn.BatchNorm2d(num_features=features),
            nn.ReLU(inplace=False),
            nn.Dropout(p=0.5)
        )

if __name__ == "__main__":
    # device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    device = "cpu"
    model = GeneratorUNet().to(device)

    input = torch.randn(2, 512, 512).to(device)
    output = model(input)
    print(output.shape)

    from torchsummary import summary

    summary(model, input_size=(512, 512))