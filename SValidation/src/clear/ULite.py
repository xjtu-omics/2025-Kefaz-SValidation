import torch
import torch.nn as nn
import torch.nn.functional as F

class AxialDW(nn.Module):
    def __init__(self, dim, mixer_kernel, dilation=1):
        super().__init__()
        h, w = mixer_kernel
        self.dw_h = nn.Conv2d(dim, dim, kernel_size=(h, 1), padding='same', groups=dim, dilation=dilation)
        self.dw_w = nn.Conv2d(dim, dim, kernel_size=(1, w), padding='same', groups=dim, dilation=dilation)



    def forward(self, x):
        # x = x + self.dw_h(x) + self.dw_w(x)
        x1 = self.dw_h(x) + self.dw_w(x)
        # X2 = self.dw_h(x) + self.dw_w(x)
        return x1

class Conv2dLayer(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size=1):
        super(Conv2dLayer, self).__init__()
        self.conv = nn.Conv2d(in_channels, out_channels, kernel_size, stride=1, padding=kernel_size//2)

    def forward(self, x):
        # x shape: (batch_size, in_channels, height, width)
        x = self.conv(x) # (batch_size, out_channels, height, width)
        return x

class SelfAttention2D(nn.Module):
    def __init__(self, in_channels, out_channels, heads=8):
        super(SelfAttention2D, self).__init__()
        self.heads = heads
        self.scale = (in_channels // heads) ** -0.5

        self.query_conv = Conv2dLayer(in_channels, out_channels)
        self.key_conv = Conv2dLayer(in_channels, out_channels)
        self.value_conv = Conv2dLayer(in_channels, out_channels)

        self.out_conv = nn.Conv2d(out_channels, out_channels, 1)

    def forward(self, x):
        batch_size, channels, height, width = x.shape
        # Apply convolutional layers
        Q = self.query_conv(x).view(batch_size, self.heads, -1, height*width)
        K = self.key_conv(x).view(batch_size, self.heads, -1, height*width)
        V = self.value_conv(x).view(batch_size, self.heads, -1, height*width)

        # Scaled Dot-Product Attention
        scores = torch.einsum('bnqd,bnkd->bnqk', Q, K) * self.scale
        attn = F.softmax(scores, dim=-1)
        context = torch.einsum('bnqk,bnkd->bnqd', attn, V).reshape(batch_size, -1, height, width)

        # Final convolutional layer
        out = self.out_conv(context)

        return out

class EncoderBlock(nn.Module):
    """Encoding then downsampling"""

    def __init__(self, in_c, out_c, mixer_kernel=(7, 7)):
        super().__init__()
        self.dw = AxialDW(in_c, mixer_kernel=(7, 7))
        self.SA = SelfAttention2D(in_c, in_c)

        self.bn = nn.BatchNorm2d(in_c)
        self.pw = nn.Conv2d(in_c, out_c, kernel_size=1)
        self.down = nn.MaxPool2d((2, 2))
        self.act = nn.GELU()

    def forward(self, x):
        # print(x.shape)
        skip = self.bn(self.dw(x))
        # skip = self.bn(self.dw(x))
        # skip = self.SA(self.dw(x))
        skip0 = self.SA(x)
        x = self.act(self.down(self.pw(skip)))
        return x, skip0


class DecoderBlock(nn.Module):
    """Upsampling then decoding"""

    def __init__(self, in_c, out_c, mixer_kernel=(7, 7)):
        super().__init__()
        self.up = nn.Upsample(scale_factor=2, mode='bilinear', align_corners=False)
        # self.up = nn.ConvTranspose2d(in_c, in_c, kernel_size=2, stride=2)
        self.pw = nn.Conv2d(in_c + out_c, out_c, kernel_size=1)
        self.bn = nn.BatchNorm2d(out_c)
        self.dw = AxialDW(out_c, mixer_kernel=(7, 7))
        self.act = nn.GELU()
        self.pw2 = nn.Conv2d(out_c, out_c, kernel_size=1)

    def forward(self, x, skip):
        x = self.up(x)
        x = torch.cat([x, skip], dim=1)
        # x = self.act(self.pw2(self.dw(self.bn(self.pw(x)))))
        x = self.act(self.pw2(self.dw(self.pw(x))))
        return x

class BottleNeckBlock(nn.Module):
    """Axial dilated DW convolution"""

    def __init__(self, dim):
        super().__init__()

        gc = dim // 4
        self.pw1 = nn.Conv2d(dim, gc, kernel_size=1)
        self.dw1 = AxialDW(gc, mixer_kernel=(3, 3), dilation=1)
        self.dw2 = AxialDW(gc, mixer_kernel=(3, 3), dilation=2)
        self.dw3 = AxialDW(gc, mixer_kernel=(3, 3), dilation=3)


        self.bn = nn.BatchNorm2d(4 * gc)
        self.pw2 = nn.Conv2d(4 * gc, dim, kernel_size=1)
        self.act = nn.GELU()

    def forward(self, x):
        x = self.pw1(x)
        x1 = self.dw1(x)
        x2 = self.dw2(x)
        x3 = self.dw3(x)
        x = torch.cat([x, x1, x2, x3], 1)
        # x = self.act(self.pw2(self.bn(x)))
        # x = self.ca(x)
        x = self.act(self.pw2(x))
        return x


class ULite(nn.Module):
    def __init__(self):
        super().__init__()

        """Encoder"""
        self.conv_in = nn.Conv2d(3, 16, kernel_size=7, padding='same')
        self.e1 = EncoderBlock(16, 32)
        self.e2 = EncoderBlock(32, 64)
        self.e3 = EncoderBlock(64, 128)
        self.e4 = EncoderBlock(128, 256)
        self.e5 = EncoderBlock(256, 512)

        """Bottle Neck"""
        self.b5 = BottleNeckBlock(512)
        # self.b5 = BottleNeckBlock(256)

        """Decoder"""
        self.d5 = DecoderBlock(512, 256)
        self.d4 = DecoderBlock(256, 128)
        self.d3 = DecoderBlock(128, 64)
        self.d2 = DecoderBlock(64, 32)
        self.d1 = DecoderBlock(32, 16)
        self.conv_out = nn.Conv2d(16, 3, kernel_size=1)

    def forward(self, x):
        # print(x.shape)
        """Encoder"""
        x = self.conv_in(x.unsqueeze(dim=0))
        x, skip1 = self.e1(x)
        x, skip2 = self.e2(x)
        x, skip3 = self.e3(x)
        x, skip4 = self.e4(x)
        x, skip5 = self.e5(x)

        """BottleNeck"""
        x = self.b5(x)  # (512, 8, 8)

        """Decoder"""
        x = self.d5(x, skip5)
        x = self.d4(x, skip4)
        x = self.d3(x, skip3)
        x = self.d2(x, skip2)
        x = self.d1(x, skip1)
        x = self.conv_out(x)
        return x


if __name__ == "__main__":
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(device)
    model = ULite().to(device)
    input = torch.randn(2, 3, 256, 256).to(device)
    output = model(input)
    print(output.shape)
    from torchsummary import summary

    summary(model, input_size=(3, 256, 256))
