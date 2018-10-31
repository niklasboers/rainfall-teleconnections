# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from pylab import *
from matplotlib.font_manager import FontProperties

def basemap(dat, lat, lon, lat_min = None, lat_max = None, lon_min = None, lon_max = None, res = 'c', proj = 'mill', origin = 0, shift = 0, title = "", loc = None, loc1 = None, loc2 = None, lsmask = False, ms_lat = None, ms_lon = None, lms = 0.1, ms_length = None, uwinds = None, vwinds = None, hatch = None, hatch_color = 'none', gc_lat0 = None, gc_lon0 = None, gc_lat1 = None, gc_lon1 = None, gc_lat2 = None, gc_lon2 = None, gdists = None, gc_lw = .5, gc_color = plt.cm.jet, gc_contours = None, gc_extend = 'both', gc_cbar = False, font_size = 8,asp = None, quiverscale = 100, stream = None, dat2 = None, dat3 = None, dat4 = None, dat5 = None, dat6 = None, area1 = None, area2 = None, area3 = None, area4 = None, area5 = None, area6 = None, area7 = None, area8 = None, area9 = None, area10 = None, area11 = None, area12 = None, area13 = None, contours = None, contours2 = None, contours3 = None, contours4 = None, contours5 = None, contours6 = None, lcontours = None, lcd = None, lclabel = None, clim = None, over = None, extend = 'neither', shaderelief = False, etopo = False, colors = None, colors2 = None, colors3 = None, colors4 = None, colors5 = None, colors6 = None,  color = plt.cm.gist_stern_r, colorbar = True, cbar_title = "",  hillshade = False, scale = .1, alpha = 1, line_color1 = None, line_color2 = None, line_color3 = None, line_color4 = None, line_color5 = None, line_color6 = None, line_color7 = None, line_color8 = None, line_color9 = None, line_color10 = None, line_color11 = None, line_color12 = None, line_color13 = None, line_style1 = 'solid', line_style2 = 'solid', line_style3 = 'solid', line_style4 = 'solid', line_style5 = 'solid', line_style6 = 'solid', line_style7 = 'solid', line_style8 = 'solid', line_style9 = 'solid', line_style10 = 'solid', line_style11 = 'solid', line_style12 = 'solid', line_style13 = 'solid', parallels = np.arange(-90., 90., 20.), meridians = np.arange(0., 360., 20.), plabels = [1, 0, 0, 0], mlabels = [0, 0, 0, 1]):
	"""plots a pcolormesh of 2-dim. dat on a geographical grid
		corresponding to lon and lat, optionally takes fixed colorlimits climl and climr, and saves the figure as title"""
	if proj == 'mill' or proj == 'cea' or proj == 'cyl':
		if lat_min is None:
			m=Basemap(projection = proj, lat_ts=12,llcrnrlon=lon.min() + shift, urcrnrlon=lon.max() + shift, llcrnrlat=lat.min(), urcrnrlat=lat.max(),resolution = res)
		else:
			m=Basemap(projection = proj, lat_ts=12,llcrnrlon=lon_min, urcrnrlon=lon_max, llcrnrlat=lat_min, urcrnrlat=lat_max,resolution = res)
		# m=Basemap(projection = proj, lat_ts=None,llcrnrlon=10, urcrnrlon=190, llcrnrlat=lat.min(), urcrnrlat=lat.max(),resolution = res)
		# m=Basemap(projection=proj, lat_0 = 0, lon_0 = shift,resolution = 'c')
	elif proj == 'ortho':
		m=Basemap(projection='ortho', lat_0 = 0., lon_0 = 0.,resolution = 'c')
	elif proj == 'moll' or proj == 'hammer' or proj == 'robin':
		m=Basemap(projection=proj, lon_0 = shift,resolution = 'c')
	elif proj == 'npstere':
		m = Basemap(projection = proj, lat_0 = 0., lon_0 = 0., boundinglat = 40, resolution = 'c')
	# dat, lons = shiftgrid(0., dat, lon)
	# print lon
	lons, dat = m.shiftdata(lon, dat, lon_0 = origin + shift)
	if lcd is not None:
		lons, lcd = m.shiftdata(lon, lcd, lon_0 = origin + shift)
	# lons = lon
	# print lons
	if loc is not None:
		xloc, yloc = m(loc[0], loc[1])
		m.plot(xloc, yloc, marker = 'D', color = 'r', markersize = lms)
	if loc1 is not None:
		xloc1, yloc1 = m(loc1[0], loc1[1])
		m.plot(xloc1, yloc1, marker = 'D', color = 'royalblue')
	if loc2 is not None:
		xloc2, yloc2 = m(loc2[0], loc2[1])
		m.plot(xloc2, yloc2, marker = 'D', color = 'darkorange')
	# dat, lons = addcyclic(dat, lons)
	if hatch is not None:
		lons, hatch = m.shiftdata(lon, hatch, lon_0 = origin + shift)
	if area1 is not None:
		lons, area1 = m.shiftdata(lon, area1, lon_0 = origin + shift)
	if area2 is not None:
		lons, area2 = m.shiftdata(lon, area2, lon_0 = origin + shift)
	if area3 is not None:
		lons, area3 = m.shiftdata(lon, area3, lon_0 = origin + shift)
	if area4 is not None:
		lons, area4 = m.shiftdata(lon, area4, lon_0 = origin + shift)
	if area5 is not None:
		lons, area5 = m.shiftdata(lon, area5, lon_0 = origin + shift)
	if area6 is not None:
		lons, area6 = m.shiftdata(lon, area6, lon_0 = origin + shift)
	if area7 is not None:
		lons, area7 = m.shiftdata(lon, area7, lon_0 = origin + shift)
	if area8 is not None:
		lons, area8 = m.shiftdata(lon, area8, lon_0 = origin + shift)
	if area9 is not None:
		lons, area9 = m.shiftdata(lon, area9, lon_0 = origin + shift)
	if area10 is not None:
		lons, area10 = m.shiftdata(lon, area10, lon_0 = origin + shift)
	if area11 is not None:
		lons, area11 = m.shiftdata(lon, area11, lon_0 = origin + shift)
	if area12 is not None:
		lons, area12 = m.shiftdata(lon, area12, lon_0 = origin + shift)
	if area13 is not None:
		lons, area13 = m.shiftdata(lon, area13, lon_0 = origin + shift)
	# m=Basemap(projection='mill', lat_0 = 0, lon_0 = 160, resolution = 'c')
	Lon, Lat = np.meshgrid(lons,lat)
	x, y = m(Lon, Lat)
	# m.drawcoastlines(linewidth = .3)
	if lsmask is True:
		m.drawlsmask(land_color = '.1', ocean_color = '1.', resolution = 'h')
	else:
		m.drawcoastlines(linewidth = .3)
		m.drawcountries(linewidth = .3)
		m.drawmapboundary()
		# m.readshapefile('/Users/omarlittle/Desktop/Desktop/PIK/climatenetwork/basemap/10m_physical/ne_10m_coastline', 'scale rank', drawbounds=True)
	# m.drawmapbound)
	# m.drawcountries(linewidth = .3)
	m.drawparallels(parallels, labels = plabels, rotation = '90', linewidth = .3, fontsize = font_size, family = 'times new roman')
	m.drawmeridians(meridians, labels = mlabels, linewidth = .3, fontsize = font_size, family = 'times new roman')
	if etopo is True:
		m.etopo()
	if shaderelief is True:
		m.shadedrelief()
	if hillshade is True:
		m.warpimage(image = '/Users/omarlittle/Desktop/Desktop/PIK/climatenetwork/basemap/SR_50M.tif', scale = scale)
	if contours is None:
		cs = m.pcolormesh(x, y, dat,  cmap = color, alpha = alpha)
		if clim is not None:
			cs.set_clim(clim)
	else:
		if colors is None:
			cs = m.contourf(x, y, dat, contours, alpha = alpha,  cmap = color, extend = extend)
		else:
			cs = m.contourf(x, y, dat, contours, alpha = alpha,  colors = colors, extend = extend)
	if contours2 is not None:
		csa = m.contourf(x, y, dat2, contours2, alpha = alpha,  colors = colors2, extend = extend)
	if contours3 is not None:
		csb = m.contourf(x, y, dat3, contours3, alpha = alpha,  colors = colors3, extend = extend)
	if contours4 is not None:
		csc = m.contourf(x, y, dat4, contours4, alpha = alpha,  colors = colors4, extend = extend)
	if contours5 is not None:
		csd = m.contourf(x, y, dat5, contours5, alpha = alpha,  colors = colors5, extend = extend)
	if contours6 is not None:
		cse = m.contourf(x, y, dat6, contours6, alpha = alpha,  colors = colors6, extend = extend)
	if hatch is not None:
		# hat = m.contourf(x, y, hatch, [-2,0, 2], colors = 'none', hatches = [None, '////'])
		mpl.rcParams['hatch.linewidth'] = 0.1
		hat = m.contourf(x, y, hatch, [-2,0, 2], colors = 'none', hatches = [None, '///////'])
	if lcontours is not None:
		lc = m.contour(x, y, lcd, lcontours, colors = 'black', linewidths = .2)
		if lclabel is True:
			plt.clabel(lc, inline = 1, fontsize = 8, fmt = '%d')
	if area1 is not None:
		cs1 = m.contour(x, y, area1, [0], linestyles = line_style1, linewidths = 1., colors = line_color1)
	if area2 is not None:
		cs2 = m.contour(x, y, area2, [0], linestyles = line_style2, linewidths = 1., colors = line_color2)
	if area3 is not None:
		cs3 = m.contour(x, y, area3, [0], linestyles = line_style3, linewidths = 1., colors = line_color3)
	if area4 is not None:
		cs4 = m.contour(x, y, area4, [0], linestyles = line_style4, linewidths = 1., colors = line_color4)
	if area5 is not None:
		cs5 = m.contour(x, y, area5, [0], linestyles = line_style5, linewidths = 1., colors = line_color5)
	if area6 is not None:
		cs6 = m.contour(x, y, area6, [0], linestyles = line_style6, linewidths = 1., colors = line_color6)
	if area7 is not None:
 		cs7 = m.contour(x, y, area7, [0], linestyles = line_style7, linewidths = 1., colors = line_color7)
	if area8 is not None:
		cs8 = m.contour(x, y, area8, [0], linestyles = line_style8, linewidths = 1., colors = line_color8)
	if area9 is not None:
		cs9 = m.contour(x, y, area9, [0], linestyles = line_style9, linewidths = 1., colors = line_color9)
	if area10 is not None:
		cs10 = m.contour(x, y, area10, [0], linestyles = line_style10, linewidths = .5, colors = line_color10)
	if area11 is not None:
		cs11 = m.contour(x, y, area11, [0], linestyles = line_style11, linewidths = .5, colors = line_color11)
	if area12 is not None:
		cs12 = m.contour(x, y, area12, [0], linestyles = line_style12, linewidths = .5, colors = line_color12)
	if area13 is not None:
		cs13 = m.contour(x, y, area13, [0], linestyles = line_style13, linewidths = .5, colors = line_color13)

        if stream is not None:
		print "be careful ..."
		lons, stream[0] = m.shiftdata(lon, stream[0], lon_0 = origin + shift)
		lons, stream[1] = m.shiftdata(lon, stream[1], lon_0 = origin + shift)
		lons, stream[2] = m.shiftdata(lon, stream[2], lon_0 = origin + shift)
		xxs, yys = m.makegrid(stream[0].shape[1], stream[0].shape[0], returnxy=True)[2:4]
		st = m.streamplot(xxs, yys, stream[1], stream[2], color="w", cmap=plt.cm.spectral, arrowsize=.5, density=3, linewidth = .5)
	FP = FontProperties()
	FP.set_family('times new roman')
	FP.set_size('%d'%font_size)
	if uwinds is not None:
		lons, uwinds = m.shiftdata(lon, uwinds, lon_0 = origin + shift)
		lons, vwinds = m.shiftdata(lon, vwinds, lon_0 = origin + shift)
		u = uwinds
		v = vwinds
		#ugrid,newlons = shiftgrid(180.,u,lon,start=False)
		#vgrid,newlons = shiftgrid(180.,v,lon,start=False)
		# transform vectors to projection grid.
		#uproj,vproj,xx,yy = m.transform_vector(ugrid,vgrid,newlons,lat,31,31,returnxy=True,masked=True)
		# now plot.
		# scale =
		Q = m.quiver(x[::asp, ::asp], y[::asp, ::asp], u[::asp, ::asp], v[::asp, ::asp],pivot = 'mid', units = 'width', scale = quiverscale)
		# make quiver key
		# qk = plt.quiverkey(Q, 0.85, 0.05, 5, '5 m/s', coordinates = 'figure', labelpos='W', labelcolor = 'black', fontproperties = {'family':'times new roman', 'size':'10', 'style':'normal'})
	if gc_lat0 is not None:
		for i in xrange(gc_lat1.shape[0]):
			if gc_color is None:
				gc = m.drawgreatcircle(gc_lon0, gc_lat0, gc_lon1[i], gc_lat1[i], linewidth = gc_lw, color = plt.cm.jet(.0001), alpha =  .02)
			else:
				gc = m.drawgreatcircle(gc_lon0, gc_lat0, gc_lon1[i], gc_lat1[i], linewidth = gc_lw, color = gc_color(int(5 * gdists[i]/gdists.max()) / 5., alpha =  .02 +  gdists[i] / gdists.max()))
	if gc_lat2 is not None:
		for i in xrange(gc_lat2.shape[0]):
			gc = m.drawgreatcircle(gc_lon0, gc_lat0, gc_lon2[i], gc_lat2[i], linewidth = gc_lw, color = plt.cm.jet(.9999), alpha =  .3)
	if loc is not None:
        	xloc, yloc = m(loc[0], loc[1])
                m.plot(xloc, yloc, marker = 'D', color = 'r', markersize = 0.5)
	if ms_lat is not None:
		m.drawmapscale(lat = ms_lat, lon = ms_lon, lon0 = 0., lat0=0., length = ms_length, barstyle='fancy')
	if colorbar is True:
		ticks_font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=font_size, weight='normal', stretch='normal')
		if gc_cbar == True:
			cs_fake = plt.contourf(x,y,np.ones_like(dat) * -1e7, gc_contours, cmap = gc_color, extend = 'neither')
			cbar = m.colorbar(cs_fake, extend = extend, orientation = 'horizontal')
		else:
			cbar = m.colorbar(cs, extend = extend)
		if cbar_title is not None:
			cbar.ax.set_ylabel(cbar_title, fontsize = font_size, family = 'times new roman')
		for label in cbar.ax.get_yticklabels():
			label.set_fontproperties(ticks_font)
	if over is not None:
		cs.cmap.set_over(over)
	plt.title('%s'%title, fontsize = font_size)
	return cs
