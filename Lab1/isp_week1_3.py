from numpy import *
import sys
import math
import cv2
import matplotlib.pyplot as plt
img1= cv2.imread("lena_translate.png",0)
img2= cv2.imread("pisa_rotate.png",0)
img3= cv2.imread("cells_scale.png",0)

a1=0.8
#a=1.3
width3, height3= img3.shape
print(width3, height3)

img_padded3 = zeros((width3+2, height3+2))
img_padded3[1:-1, 1:-1] = img3
img3_s=zeros((width3, height3))

xc=width3/2
yc=height3/2
for xt in range(width3):
  for yt in range(height3):
    xs = (xt-xc)/a1 + xc
    ys = (yt-yc)/a1+yc 
    x=xs+1
    y=ys+1
    xs1=math.floor(x)
    ys1=math.floor(y)
    a=x-xs1
    b=y-ys1
    if xs1>=0 and xs1<=width3 and ys1>=0 and ys1<=height3:
      img3_s[xt,yt]=(1-a)*(1-b)*img_padded3[xs1,ys1]+(1-a)*b*img_padded3[xs1,ys1+1]+a*(1-b)*img_padded3[xs1+1,ys1]+a*b*img_padded3[xs1+1,ys1+1]
    else:
      img3_s[xt,yt]=255

#fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,12), constrained_layout=True)
#ax1.imshow(img3,'gray')
#ax1.title.set_text("Actual Image")
#ax2.imshow(img3_s,'gray')
#ax2.title.set_text("Scaled Image by a factor 0.8")
