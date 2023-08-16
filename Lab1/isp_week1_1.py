from numpy import *
import sys
import math
import cv2
import matplotlib.pyplot as plt
img1= cv2.imread("lena_translate.png",0)
img2= cv2.imread("pisa_rotate.png",0)
img3= cv2.imread("cells_scale.png",0)


tx=3.75
ty=4.3
width, height= img1.shape
#print(width, height)
img1_t=zeros((width, height))
img_padded = zeros((width+2, height+2))
img_padded[1:-1, 1:-1] = img1 

for xt in range(width):
  for yt in range(height):
    xs=xt-tx
    ys=yt-ty
    #intensity=img1[xs,ys]
    x=xs+1
    y=ys+1
    xs1=math.floor(x)
    ys1=math.floor(y)
    a=x-xs1
    b=y-ys1
    if xs1>=0 and xs1<=width and ys1>=0 and ys1<=height:
      img1_t[xt,yt]=(1-a)*(1-b)*img_padded[xs1,ys1]+(1-a)*b*img_padded[xs1,ys1+1]+a*(1-b)*img_padded[xs1+1,ys1]+a*b*img_padded[xs1+1,ys1+1]
    else:
      img1_t[xt,yt]=0

#plt.figure(1)
#plt.imshow(img1,'gray')
#savefig("222.png")
#plt.show()
cv2.imshow('image', img1)
savefig('kl.png')

#fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,12), constrained_layout=True)
#ax1.imshow(img1,'gray')
#ax2.imshow(img1_t,'gray')
#plt.show()
#print(img1_t)
