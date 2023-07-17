# 包中image项目中的展示图像的几个函数
# log对比度是按照ds9的log公式来的，不太一样
import numpy as np

def logstretch(img,a=None):
	# img: 数据，二维矩阵形式
	# a: log对比度的参数，光学图像典型值一般用1000，
	# 	红外典型值一般用100，越高能的波段用的值应该越高，100和1000只是参考值
	if a==None:
		a = 500
	img_log = np.log10(a * (img - np.min(img)) + 1) / np.log10(a) #img - np.min(img)和+1是为了不要造成log的时候成负数
	return img_log

def rescale(img):
	# img: 数据，二维矩阵形式
	img_rescale = (img - np.amin(img)) / (np.amax(img) - np.amin(img))
	return img_rescale

def positivy(img):
	# img: 数据，二维矩阵形式
	img_rescale = 0.001 * (np.amax(img) - np.amin(img)) + img - np.amin(img)
	return img_rescale