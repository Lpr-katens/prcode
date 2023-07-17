# ===================================================
# 此函数用来平滑光谱，可以用均值或者中值
# ===================================================
import numpy as np

def average(flux,pixels=None):
	# ------------------------------------------------
	# 输入：
	# 	流量，以及想做平均的前后pixels个像素的值
	# ------------------------------------------------
	if pixels == None:
		pixs = 3
	else:
		pixs = pixels
	# 生成pixs个 flux长度+pixs-1 的数组，即采用数组前挪(pixs-1)/2与后挪(pixs-1)/2个单位，然后pixs个数组的每个位置对应叠加取均值的方法
	array = np.full((pixs,flux.size + pixs-1),0.)
	# 先把原flux放在矩阵的第一行
	array[0,int((pixs-1) / 2):-int((pixs-1) / 2)] = flux
	# 分别向前挪动flux 1……(pixs-1)/2 个单位，和向后挪动 1……(pixs-1)/2 个单位
	for num in range(1,int((pixs-1) / 2) + 1):
		array[num,int((pixs-1) / 2) - num:-int((pixs-1) / 2) - num] = flux
		array[num + int((pixs-1) / 2),int((pixs-1) / 2) + num:int((pixs-1) / 2) + num + flux.size] = flux
	# 求平均，其实就是对每个pixel求前后(pixs-1)/2个pixels的均值，来替换这个pixel的值
	flux_aver = np.mean(array,axis = 0)
	return flux_aver[int((pixs-1) / 2):-int((pixs-1) / 2)]

def median(flux,pixels=None):
	# ------------------------------------------------
	# 输入：
	# 	流量，以及想求中值的前后pixels个像素的值
	# 跟函数average一样的思路，只不过最后取的是中值而非均值
	# ------------------------------------------------
	if pixels == None:
		pixs = 3
	else:
		pixs = pixels
	# 生成pixs个 flux长度+pixs-1 的数组，即采用数组前挪(pixs-1)/2与后挪(pixs-1)/2个单位，然后pixs个数组的每个位置对应叠加取均值的方法
	array = np.full((pixs,flux.size + pixs-1),0.)
	# 先把原flux放在矩阵的第一行
	array[0,int((pixs-1) / 2):-int((pixs-1) / 2)] = flux
	# 分别向前挪动flux 1……(pixs-1)/2 个单位，和向后挪动 1……(pixs-1)/2 个单位
	for num in range(1,int((pixs-1) / 2) + 1):
		array[num,int((pixs-1) / 2) - num:-int((pixs-1) / 2) - num] = flux
		array[num + int((pixs-1) / 2),int((pixs-1) / 2) + num:int((pixs-1) / 2) + num + flux.size] = flux
	# 求中值，其实就是对每个pixel求前后(pixs-1)/2个pixels的中值，来替换这个pixel的值
	flux_aver = np.median(array,axis = 0)
	return flux_aver[int((pixs-1) / 2):-int((pixs-1) / 2)]