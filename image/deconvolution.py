import numpy as np

def deconvolved_color(array, kernel):
    # array: 观测图像
    # kernel: 红端到蓝端观测图像的psf的kernel

    if array.shape[0] < kernel.shape[0] or array.shape[1] < kernel.shape[1]:
        raise ValueError('Kernel must have smaller size than image!')
    
    image = array.copy()
    kernel = kernel.copy()
    # 定义一个矩阵matrix，用来解方程，大小为(image.shape[0] * image.shape[1], image.shape[0] * image.shape[1])，也就是图像总像素数的平方
    matrix = np.zeros((image.shape[0] * image.shape[1], image.shape[0] * image.shape[1]))
    kernel_x_size, kernel_y_size = [int((kernel.shape[0] - 1) / 2), int((kernel.shape[1] - 1) / 2)]
    img_x_size, img_y_size = [image.shape[0], image.shape[1]]
    
    # 开始为矩阵matrix的每一行填数，方程的右侧列向量为图像image从 第一行0:结尾->第二行0:结尾 一直到最后一个数，自变量的列向量为待求的卷积psf之前的高清图像 从第一行0:结尾->第二行0:结尾 一直到最后一个数
    for i in range(img_x_size):
        for j in range(img_y_size):
            # 每次循环都要重新定义matrix_image，用来计算matrix的第i×img_y_size+j行应该是什么样
            matrix_image = np.zeros((img_x_size, img_y_size))
            # 把第i,j个像素卷积的时候所涉及到的哪些像素都赋值上kernel的值，其余的值在卷积的时候都与第i,j个像素无关，所以都设置为0
            if i - kernel_x_size < 0 or j - kernel_y_size < 0:
                matrix_image[max(0, i - kernel_x_size) : min(i + kernel_x_size + 1, img_x_size), max(0, j - kernel_y_size) : min(j + kernel_y_size + 1, img_y_size)] = kernel[{True: abs(i - kernel_x_size), False: 0}[i - kernel_x_size < 0] : {True: -(i + kernel_x_size + 1 - img_x_size), False: kernel_x_size * 2 + 1}[i + kernel_x_size > img_x_size - 1], {True: abs(j - kernel_y_size), False: 0}[j - kernel_y_size < 0] : {True: -(j + kernel_y_size + 1 - img_y_size), False: kernel_y_size * 2 + 1}[j + kernel_y_size > img_y_size - 1]]
            elif i + kernel_x_size > img_x_size - 1 or j + kernel_y_size > img_y_size - 1:
                matrix_image[i - kernel_x_size : min(i + kernel_x_size + 1, img_x_size), j - kernel_y_size : min(j + kernel_y_size + 1, img_y_size)] = kernel[0 : {True: min(-(kernel_x_size - (img_x_size - i - 1)), -1), False: kernel_x_size * 2 + 1}[i + kernel_x_size > img_x_size - 1], 0 : {True: min(-(kernel_y_size - (img_y_size - j - 1)), -1), False: kernel_y_size * 2 + 1}[j + kernel_y_size > img_y_size - 1]]
            else:
                matrix_image[i - kernel_x_size : i + kernel_x_size + 1, j - kernel_y_size : j + kernel_y_size + 1] = kernel
            # 把这个第i,j个像素在矩阵matrix对应的i×img_y_size+j行，赋上值为matrix_image.flatten()
            matrix[i * img_y_size + j, :] = matrix_image.flatten()
    # 矩阵matrix的逆矩阵为matrix_invert
    matrix_invert = np.linalg.inv(matrix)
    # 待求的卷积psf之前的高清图像 等于 逆矩阵matrix_invert与观测图像image.flatten()矩阵乘
    image_deconvolved = np.dot(matrix_invert, image.flatten()).reshape(img_x_size, img_y_size)
    return (image_deconvolved, matrix, matrix_invert)

# # 用cc老师的例子测试，成功
# from astropy.convolution import convolve, Gaussian2DKernel
# from astropy.modeling.models import Gaussian2D
# import matplotlib.pyplot as plt

# kernel = Gaussian2DKernel(1, x_size=3, y_size=3)
# kernel = kernel.array
# kernel = kernel/np.sum(kernel)

# xx = np.arange(10)
# yy = np.arange(10)
# x, y = np.meshgrid(xx, yy)

# fig = plt.figure(figsize=(5 * 3, 5 * 2))
# gs = fig.add_gridspec(2, 3)
# # 做个假星系
# galaxy = Gaussian2D(1, 4, 3, 2, 1, theta=0.5)
# image = galaxy(x, y)*10
# ax0 = fig.add_subplot(gs[0,0])
# img0 = ax0.imshow(image, cmap='coolwarm')
# ax0.set_title('mock galaxy')
# fig.colorbar(img0, ax=ax0)
# # 做这个假星系在有psf的影响和噪音情况下的观测图像
# conv_image = convolve(image, kernel, normalize_kernel=True)
# rms = 0.0*np.random.randn(10, 10)
# obs_image = conv_image + rms
# ax1 = fig.add_subplot(gs[0,1])
# img1 = ax1.imshow(obs_image, cmap='coolwarm')
# ax1.set_title('mock galaxy in observation')
# fig.colorbar(img1, ax=ax1)
# # 这个星系用到的 用来解方程的矩阵matrix
# matrix = deconvolved_color(obs_image, kernel)[1]
# ax2 = fig.add_subplot(gs[0,2])
# img2 = ax2.imshow(matrix, cmap='coolwarm')
# ax2.set_title('matrix')
# fig.colorbar(img2, ax=ax2)
# # 这个星系用到的 用来解方程的矩阵的逆矩阵matrix_invert
# matrix_invert = deconvolved_color(obs_image, kernel)[2]
# ax3 = fig.add_subplot(gs[1,0])
# img3 = ax3.imshow(matrix_invert, cmap='coolwarm')
# ax3.set_title('matrix invert')
# fig.colorbar(img3, ax=ax3)
# # 这个星系反卷积之后的高清图像
# image_theroy = deconvolved_color(obs_image, kernel)[0]
# ax4 = fig.add_subplot(gs[1,1])
# img4 = ax4.imshow(image_theroy, cmap='coolwarm')
# ax4.set_title('deconvolved observed mock galaxy')
# fig.colorbar(img4, ax=ax4)
# # 残差图像
# residual = image - deconvolved_color(obs_image, kernel)[0]
# ax5 = fig.add_subplot(gs[1,2])
# img5 = ax5.imshow(residual, cmap='coolwarm')
# ax5.set_title('residual(mock - deconvolved obs)')
# fig.colorbar(img5, ax=ax5)
# plt.savefig('/Users/lpr/Data/jwst_lirg_project/plots/deconvolved.png')