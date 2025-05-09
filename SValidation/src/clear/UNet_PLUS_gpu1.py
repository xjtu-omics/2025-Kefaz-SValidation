import torch
import torch.nn as nn

class encoder_conv(nn.Module):
    def __init__(self, in_chs, out_chs):
        super(encoder_conv, self).__init__()
        #self.conv1 = nn.Conv2d(in_chs, out_chs, kernel_size=3, stride=1, padding="same")
        #self.conv2 = nn.Conv2d(out_chs, out_chs, kernel_size=3, stride=1, padding="same")
        out_channels = int(out_chs / 4)
        self.maxpool = nn.MaxPool2d(kernel_size=2, stride=2)
        self.avgpool = nn.AvgPool2d(kernel_size=2, stride=2)
        self.lrelu = nn.LeakyReLU()
        #self.bn = nn.BatchNorm2d(out_chs)
        self.conv1 = nn.Conv2d(in_chs, out_channels, kernel_size=1, stride=1, dilation=1, padding="same")
        self.conv2 = nn.Conv2d(in_chs, out_channels, kernel_size=3, stride=1, dilation=2, padding="same")
        self.conv3 = nn.Conv2d(in_chs, out_channels, kernel_size=5, stride=1, dilation=3, padding="same")
        self.conv4 = nn.Conv2d(in_chs, out_channels, kernel_size=7, stride=1, dilation=3, padding="same")

    def forward(self, x):
        #x1 = self.lrelu(self.conv1(x))
        #x2 = self.lrelu(self.conv2(x1))
        #x3 = self.pool(x2)
        x1 = self.lrelu(self.conv1(x))
        x2 = self.lrelu(self.conv2(x))
        x3 = self.lrelu(self.conv3(x))
        x4 = self.lrelu(self.conv4(x))
        x_out = torch.cat([x1,x2,x3,x4], dim=1)
        out = self.maxpool(x_out) + self.avgpool(x_out)
        return out
class self_attention_gate(nn.Module):
    def __init__(self, in_chs, out_chs):
        super(self_attention_gate, self).__init__()
        self.avgpool = nn.AdaptiveAvgPool2d(1)
        self.maxpool = nn.AdaptiveMaxPool2d(1)
        self.conv_sq = nn.Conv2d(in_chs, out_chs // 16, kernel_size=1, stride=1)
        self.lrelu = nn.LeakyReLU()
        self.conv_ex = nn.Conv2d(out_chs // 16, out_chs, kernel_size=1, stride=1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, skip):
        #x1 = self.avgpool(x) + self
        skip1 = self.avgpool(skip) + self.maxpool(skip)
        w = self.sigmoid(self.conv_ex(self.lrelu(self.conv_sq(skip1))))
        return skip * w

'''class self_attention_gate_3in(nn.Module):
    def __init__(self, in_chs, out_chs):
        super(self_attention_gate_3in, self).__init__()
        self.avgpool = nn.AdaptiveAvgPool2d(1)
        self.maxpool = nn.AdaptiveMaxPool2d(1)
        self.conv_sq = nn.Conv2d(in_chs, out_chs // 16, kernel_size=1, stride=1)
        self.lrelu = nn.LeakyReLU()
        self.conv_ex = nn.Conv2d(out_chs // 16, out_chs, kernel_size=1, stride=1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, skip):
        #x1 = self.avgpool(x)
        skip1 = self.avgpool(skip) + self.maxpool(skip)
        w = self.sigmoid(self.conv_ex(self.lrelu(self.conv_sq(skip1))))
        return skip * w'''




class decoder_conv(nn.Module):
    def __init__(self, in_chs, out_chs):
        super(decoder_conv, self).__init__()
        self.SAG = self_attention_gate(out_chs, out_chs)
        self.deconv = nn.ConvTranspose2d(in_chs, out_chs, kernel_size=2, stride=2)
        self.conv1 = nn.Conv2d(in_chs, out_chs, kernel_size=3, stride=1, padding="same")
        self.conv2 = nn.Conv2d(out_chs, out_chs, kernel_size=3, stride=1, padding="same")
        self.lrelu = nn.LeakyReLU()
        self.dropout = nn.Dropout(0.3)
    def forward(self, x, skip):
        x1 = self.deconv(x)
        skip = self.SAG(skip)
        x2 = torch.cat([x1, skip], dim=1)
        x3 = self.lrelu(self.conv1(x2))
        x4 = self.lrelu(self.conv2(x3))
        return x4
class decoder_conv_3input(nn.Module):
    def __init__(self, in_chs, m_chs, out_chs):
        super(decoder_conv_3input, self).__init__()
        self.deconv = nn.ConvTranspose2d(in_chs, m_chs, kernel_size=2, stride=2)
        self.conv1 = nn.Conv2d(m_chs * 3, out_chs, kernel_size=3, stride=1, padding="same")
        self.conv2 = nn.Conv2d(out_chs, out_chs, kernel_size=3, stride=1, padding="same")
        self.lrelu = nn.LeakyReLU()
        self.SAG = self_attention_gate(out_chs, out_chs)
        #self.dropout = nn.Dropout(0.3)

    def forward(self, x, skip1, skip2):
        x1 = self.deconv(x)
        skip1 = self.SAG(skip1)
        skip2 = self.SAG(skip2)
        x2 = torch.cat([x1, skip1, skip2], dim=1)
        x3 = self.lrelu(self.conv1(x2))
        x4 = self.lrelu(self.conv2(x3))
        return x4

'''class ASPPBlock(nn.Module):
    def __init__(self, in_channels, out_channels_list, kernel_size_list, dilation_list):
        super(ASPPBlock, self).__init__()
        self.conv_num = len(out_channels_list)
        assert (self.conv_num == 4)
        assert (self.conv_num == len(kernel_size_list) and self.conv_num == len(dilation_list))

        self.conv_1 = nn.Conv2d(in_channels, out_channels_list[0], kernel_size=kernel_size_list[0],
                                dilation=dilation_list[0], padding='same')
        self.conv_2 = nn.Conv2d(in_channels, out_channels_list[1], kernel_size=kernel_size_list[1],
                                dilation=dilation_list[1], padding='same')
        self.conv_3 = nn.Conv2d(in_channels, out_channels_list[2], kernel_size=kernel_size_list[2],
                                dilation=dilation_list[2], padding='same')
        self.conv_4 = nn.Conv2d(in_channels, out_channels_list[3], kernel_size=kernel_size_list[3],
                                dilation=dilation_list[3], padding='same')

        out_channels = out_channels_list[0] + out_channels_list[1] + out_channels_list[2] + out_channels_list[3]
        self.conv_1x1 = nn.Sequential(
            nn.Conv2d(out_channels, out_channels, kernel_size=1, padding=0),
            nn.BatchNorm2d(out_channels),
            nn.LeakyReLU())

    def forward(self, x):
        x1 = self.conv_1(x)
        x2 = self.conv_2(x)
        x3 = self.conv_3(x)
        x4 = self.conv_4(x)

        y = torch.cat([x1, x2, x3, x4], dim=1)
        y = self.conv_1x1(y)
        return y'''



class UNet_PP(nn.Module):
    def __init__(self, in_chs=3):
        super(UNet_PP,self).__init__()
        self.encoder1 = encoder_conv(in_chs, 64)
        self.encoder2 = encoder_conv(64, 128)
        self.encoder3 = encoder_conv(128, 256)
        self.encoder4 = encoder_conv(256, 512)
        self.decoder1 = decoder_conv(128, 64)
        self.decoder2 = decoder_conv(256, 128)
        self.decoder3 = decoder_conv(512, 256)
        self.decoder4 = decoder_conv_3input(128, 64, 64)
        self.decoder5 = decoder_conv_3input(256, 128, 128)
        self.decoder6 = decoder_conv_3input(128, 64, 64)
        self.outconv = nn.ConvTranspose2d(64 * 3, 3, kernel_size=2, stride=2)

        '''f4 = 512
        aspp_chns = [int(f4 / 4), int(f4 / 4), int(f4 / 4), int(f4 / 4)]
        aspp_knls = [1, 3, 3, 3]
        aspp_dila = [1, 2, 4, 6]
        self.aspp = ASPPBlock(f4, aspp_chns, aspp_knls, aspp_dila)'''

    def forward(self, x):
        x = x.unsqueeze(dim=0)
        x1_0 = self.encoder1(x)
        x2_0 = self.encoder2(x1_0)
        x3_0 = self.encoder3(x2_0)
        x4_0 = self.encoder4(x3_0)
        #x4_0 = self.aspp(x4_0)
        x1_1 = self.decoder1(x2_0, x1_0)
        x2_1 = self.decoder2(x3_0, x2_0)
        x3_1 = self.decoder3(x4_0, x3_0)

        x1_2 = self.decoder4(x2_1, x1_1, x1_0)
        x2_2 = self.decoder5(x3_1, x2_1, x2_0)

        x1_3 = self.decoder6(x2_2, x1_1, x1_0)
        x_out = torch.cat([x1_1, x1_2, x1_3], dim=1)
        out = self.outconv(x_out)
        return out

if __name__ == "__main__":
    from torchsummary import summary
    model = UNet_PP().to('cuda')
    input = torch.randn(1, 3, 256, 256).to('cuda')
    out = model(input)
    summary(model, input_size=(3, 256, 256))
    print(out.shape)

