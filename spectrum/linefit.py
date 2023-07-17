# ===================================================
# 此函数用来拟合发射线或吸收线，包括连续谱的拟合
# ===================================================
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats

# 定义高斯函数
def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(- (x-x0) ** 2 / (2 * sigma**2))

# 定义带斜度的连续谱+高斯函数
def gaussPstraightline(x,H,A,x0,sigma,k):
	gaussian = H + A*np.exp(- (x-x0) ** 2 / (2 * sigma**2))
	baseline = k*x
	return gaussian+baseline

# 定义双高斯函数
def Ha(x,H,A1,x0,sigma1,A2,sigma2):
	gaussian1 = A1*np.exp(- (x-x0) ** 2 / (2 * sigma1**2))
	gaussian2 = A2*np.exp(- (x-x0) ** 2 / (2 * sigma2**2))
	f = gaussian1+gaussian2+H
	return f

# 定义四高斯函数
def HaPNII(x,H,A1,x01,sigma1,A2,sigma2,A3,x03,sigma3,x04,sigma4):
	gaussian1 = A1*np.exp(- (x-x01) ** 2 / (2 * sigma1**2))
	gaussian2 = A2*np.exp(- (x-x01) ** 2 / (2 * sigma2**2))
	gaussian3 = A3*np.exp(- (x-x03) ** 2 / (2 * sigma3**2))
	# NII_ratio来自全局变量，是从HaPNII_fit来的，由用户输入
	gaussian4 = A3*NII_ratio*np.exp(- (x-x04) ** 2 / (2 * sigma4**2))
	f = gaussian1+gaussian2+gaussian3+gaussian4+H
	return f


def gauss_fit(wave,flux,peak_wave,p0=None,bounds=None):
	# ------------------------------------------------
	# 输入：
	# 	波长和对应的流量，以及在这个星系的红移下，想拟合的线的
	# 	顶峰波长是多少(就是高斯的期望值)
	# ------------------------------------------------
	if bounds == None:
		bounds = [-np.inf,-np.amax(abs(flux)),peak_wave-20,-np.inf],[np.inf,np.amax(abs(flux)),peak_wave+20,np.inf]
	else:
		bounds = bounds
	if p0 == None:
		p0 = [1,np.median(flux),peak_wave,1]
	else:
		p0 = p0
	popt,pcov = curve_fit(gauss,wave,flux,p0,bounds=(bounds))
	intercept = popt[0]
	amplitude = popt[1]
	mu = popt[2]
	standard_deviation = popt[3]
	return [intercept,amplitude,mu,standard_deviation]

def continuum_fit(wave,flux):
	# ------------------------------------------------
	# 输入：
	# 	想拟合连续谱的波长和对应的流量
	# ------------------------------------------------
	result = stats.linregress(wave,flux)
	k = result.slope
	b = result.intercept
	return [k,b]

def gaussPcontinuum_fit(wave,flux,peak_wave,p0=None,bounds=None):
	# ------------------------------------------------
	# 输入：
	# 	想拟合线+连续谱的波长和对应的流量
	# ------------------------------------------------
	if bounds == None:
		bounds = [-np.inf,-np.amax(abs(flux)),peak_wave-20,-np.inf,-np.inf],[np.inf,np.amax(abs(flux)),peak_wave+20,np.inf,np.inf]
	else:
		bounds = bounds
	if p0 == None:
		p0 = [1,np.median(flux),peak_wave,1,1.e-5]
	else:
		p0 = p0
	popt,pcov = curve_fit(gaussPstraightline,wave,flux,p0=p0,bounds=(bounds))
	intercept = popt[0]
	amplitude = popt[1]
	mu = popt[2]
	standard_deviation = popt[3]
	k = popt[4]
	return [intercept,amplitude,mu,standard_deviation,k]

def Ha_fit(wave,flux,peak_wave,p0=None,bounds=None):
	# ----------------------------------------------
	# 输入：
	# 	光谱的波长和流量，ha的波长，窄线+宽线
	# 	1.宽线的展宽不小于500km/s
	# 	2.必须换算成静止系波长再做拟合
	# 	3.必须减掉连续谱之后再拟合
	# ----------------------------------------------
	sigma_broad_lolimit = 5.*peak_wave / (3.e3* 2*np.sqrt(2*np.log(2)))
	peak_flux = np.amax(abs(flux[(wave>peak_wave-40)&(wave<peak_wave+40)]))
	if bounds == None:
		bounds = [-np.inf,-peak_flux,peak_wave-20,0,-peak_flux,sigma_broad_lolimit],[np.inf,peak_flux,peak_wave+20,np.inf,peak_flux,np.inf]
	else:
		bounds = bounds
	
	if p0 == None:
		p0 = [0,peak_flux,peak_wave,30,peak_flux,peak_wave,sigma_broad_lolimit]
	else:
		p0 = p0
	popt,pcov = curve_fit(Ha,wave,flux,p0=p0,bounds=(bounds))
	intercept = popt[0]
	amplitude_narrow = popt[1]
	mu_narrow = popt[2]
	standard_deviation_narrow = popt[3]
	amplitude_broad = popt[4]
	standard_deviation_broad = popt[5]
	return [intercept,amplitude_narrow,mu_narrow,standard_deviation_narrow,amplitude_broad,standard_deviation_broad]

def HaPNII_fit(wave,flux,peak_wave_Ha,peak_wave_NII,NII_flux_ratio=None,p0=None,bounds=None):
	# ----------------------------------------------
	# 输入：
	# 	光谱的波长和流量，ha的波长，窄线+宽线
	# 	1.宽线的展宽不小于500km/s
	# 	2.必须换算成静止系波长再做拟合
	# 	3.必须减掉连续谱之后再拟合
	# 	4.如果NII线比没有规定，那么就用Greene+05中的值2.96
	# ----------------------------------------------
	sigma_broad_lolimit = 5.*peak_wave_Ha / (3.e3* 2*np.sqrt(2*np.log(2))) #宽线最窄500km/s
	sigma_broad_uplimit = 100.*peak_wave_Ha / (3.e3* 2*np.sqrt(2*np.log(2))) #宽线最宽10000km/s
	peak_flux = np.amax(abs(flux[(wave>min(peak_wave_NII)-40)&(wave<max(peak_wave_NII)+40)]))
	if NII_flux_ratio == None:
		NII_flux_ratio = 2.96 # 2.96 from Greene+05
	# 定义NII的全局变量
	global NII_ratio
	NII_ratio = NII_flux_ratio
	# peak_wave_NII需要包含短波和长波两个中心波长
	if len(peak_wave_NII) < 2:
		raise ValueError('NII must have two central wavelength')
	if bounds == None:
		bounds = [-np.inf,0,peak_wave_Ha-20,0,0,sigma_broad_lolimit,0,min(peak_wave_NII)-20,0,max(peak_wave_NII)-20,0],[peak_flux,peak_flux,peak_wave_Ha+20,sigma_broad_lolimit,peak_flux,sigma_broad_uplimit,peak_flux/NII_flux_ratio,min(peak_wave_NII)+20,sigma_broad_lolimit,max(peak_wave_NII)+20,sigma_broad_lolimit]
	else:
		bounds = bounds
	if p0 == None:
		p0 = [0,peak_flux,peak_wave_Ha,1,peak_flux,sigma_broad_lolimit,peak_flux/NII_flux_ratio,min(peak_wave_NII),1,max(peak_wave_NII),1]
	else:
		p0 = p0
	popt,pcov = curve_fit(HaPNII,wave,flux,bounds=(bounds))#p0=p0,
	intercept = popt[0]
	amplitude_narrow = popt[1]
	mu_narrow = popt[2]
	standard_deviation_narrow = popt[3]
	amplitude_broad = popt[4]
	standard_deviation_broad = popt[5]
	amplitude_NIIshort = popt[6]
	mu_NIIshort = popt[7]
	standard_deviation_NIIshort = popt[8]
	mu_NIIlong = popt[9]
	standard_deviation_NIIlong = popt[10]
	return [intercept,amplitude_narrow,mu_narrow,standard_deviation_narrow,amplitude_broad,standard_deviation_broad,amplitude_NIIshort,mu_NIIshort,standard_deviation_NIIshort,mu_NIIlong,standard_deviation_NIIlong]