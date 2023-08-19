"""
Original script from Robert Barkhouser:


total_flux = 0.0
x_mean = 0.0
y_mean = 0.0
x_m2 = 0.0
y_m2 = 0.0
#
for row in range(ylb, yub, 1):
    for col in range(xlb, xub, 1):
        pxl_flux = dd[col,row]
        if pxl_flux > 0.0:
            tmp = total_flux + pxl_flux
            x = (col - dth_ctr + 0.5) * dth_scl
            y = (row - dth_ctr + 0.5) * dth_scl
            delta_x = x - x_mean
            delta_y = y - y_mean
            Rx = delta_x * pxl_flux / tmp
            Ry = delta_y * pxl_flux / tmp
            x_mean = x_mean + Rx
            y_mean = y_mean + Ry
            x_m2 = x_m2 + (total_flux * delta_x * Rx)
            y_m2 = y_m2 + (total_flux * delta_y * Ry)
            total_flux = tmp
#
xvar = x_m2 / total_flux
yvar = y_m2 / total_flux
x_sigma = sqrt(xvar)
y_sigma = sqrt(yvar)
rms = sqrt((x_sigma * x_sigma) + (y_sigma * y_sigma))

"""

import math 

total_flux = 0.0
x_mean = 0.0
y_mean = 0.0
x_m2 = 0.0
y_m2 = 0.0

for row in range(ylb, yub, 1):
    for col in range(xlb, xub, 1):
        pxl_flux = dd[row, col]
        if pxl_flux > 0.0:
            tmp = total_flux + pxl_flux
            x = col - x_ctr
            y = row - y_ctr
            delta_x = x - x_mean
            delta_y = y - y_mean
            Rx = delta_x * pxl_flux / tmp
            Ry = delta_y * pxl_flux / tmp
            x_mean = x_mean + Rx
            y_mean = y_mean + Ry
            x_m2 = x_m2 + (total_flux * delta_x * Rx)
            y_m2 = y_m2 + (total_flux * delta_y * Ry)
            total_flux = tmp

xvar = x_m2 / total_flux
yvar = y_m2 / total_flux
x_sigma = math.sqrt(xvar)
y_sigma = math.sqrt(yvar)
rms = math.sqrt((x_sigma * x_sigma) + (y_sigma * y_sigma))
