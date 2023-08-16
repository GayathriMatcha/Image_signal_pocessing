from numpy import *
import sys
import math
import cv2
import matplotlib.pyplot as plt
img1= cv2.imread("lena_translate.png",0)
img2= cv2.imread("pisa_rotate.png",0)
img3= cv2.imread("cells_scale.png",0)

#let theta=10',5',4'
theta=-5
theta_r=pi*theta/180
width2, height2= img2.shape
print(width2, height2)

img_padded2 = zeros((width2+2, height2+2))
img_padded2[1:-1, 1:-1] = img2
img2_r=zeros((width2, height2))

for xt in range(width2):
  for yt in range(height2):
      xc, yc = xt-width2/2, yt-height2/2
      xs = cos(theta_r)*xc - sin(theta_r)*yc + width2/2
      ys = cos(theta_r)*yc + sin(theta_r)*xc + height2/2
      x=xs+1
      y=ys+1
      xs1=math.floor(x)
      ys1=math.floor(y)
      a=x-xs1
      b=y-ys1
      if xs1>=0 and xs1<=width2 and ys1>=0 and ys1<=height2:
        img2_r[xt,yt]=(1-a)*(1-b)*img_padded2[xs1,ys1]+(1-a)*b*img_padded2[xs1,ys1+1]+a*(1-b)*img_padded2[xs1+1,ys1]+a*b*img_padded2[xs1+1,ys1+1]
      else:
        img2_r[xt,yt]=125

#fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,12), constrained_layout=True)
#ax1.imshow(img2,'gray')
#ax1.title.set_text("Actual Image")

#ax2.imshow(img2_r,'gray')
#ax1.title.set_text("Rotated Image") 
