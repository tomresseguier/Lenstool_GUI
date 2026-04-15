import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
from tqdm import tqdm

import PyQt5
import pyqtgraph as pg
from astropy.table import Table
import pandas as pd

###############################################################################
#from .source_extraction.match_cat import run_match
from .utils.utils_astro.cat_manip import match_cat2
from .utils.utils_plots.plot_utils_general import *
from .utils.utils_Qt.selectable_classes import SelectableEllipse, SelectableScatter, SelectSources, ellipse_maker_ROI
from .utils.utils_Qt.utils_general import make_handles, InRectangle, make_full_color
from .utils.utils_general.utils_general import make_colnames_dict






def open_cat(cat_path) :
    if '.fits' in cat_path :
        cat = Table.read(cat_path, format='fits')
    else :
        with open(cat_path, 'r') as raw_cat :
            first_line, second_line = raw_cat.readlines()[0:2]
            #print(first_line)
            start_line = 1 if second_line.startswith('--') else 0
        if len(first_line.split()) > len(first_line.split(',')) :
            cat_df = pd.read_csv(cat_path, sep='\s+', skip_blank_lines=True, comment='#')[start_line:].apply(pd.to_numeric, errors='coerce')
        else :
            cat_df = pd.read_csv(cat_path, skip_blank_lines=True, comment='#')[start_line:].apply(pd.to_numeric, errors='coerce')
        cat = Table.from_pandas(cat_df)
    return cat

    
def make_uniform_names_cat(cat, self) :
    uniform_names_cat = cat.copy()
    colnames_dict = make_colnames_dict(cat, use_default_names=self.use_default_names)
    
    self._vprint('Column names to be used:')
    self._vprint(colnames_dict)
    
    for colname in colnames_dict.keys() :
        if colnames_dict[colname] is not None :
            
            if colname in uniform_names_cat.colnames :
                uniform_names_cat[colname] = uniform_names_cat[colnames_dict[colname]]
            elif colname in ['ra', 'dec'] :
                uniform_names_cat.rename_column(colnames_dict[colname], colname)
            else :
                uniform_names_cat.add_column( uniform_names_cat[colnames_dict[colname]], name=colname )
                
            
    if colnames_dict['a'] != 'A_IMAGE' and colnames_dict['b'] != 'B_IMAGE' :
        if self.units is None :
            self.units = input("ellipticity parameters " + str(colnames_dict['a']) \
                               + " and " + str(colnames_dict['b']) + " in pixels? [y][arcsec][deg]")
        if self.units == 'deg' :
            uniform_names_cat.replace_column( 'a', uniform_names_cat['a']/(self.fits_image.pix_deg_scale) )
            uniform_names_cat.replace_column( 'b', uniform_names_cat['b']/(self.fits_image.pix_deg_scale) )
        if self.units == 'arcsec' :
            uniform_names_cat.replace_column( 'a', uniform_names_cat['a']/(self.fits_image.pix_deg_scale*3600) )
            uniform_names_cat.replace_column( 'b', uniform_names_cat['b']/(self.fits_image.pix_deg_scale*3600) )
    
    #if colnames_dict['x']==None :
    x, y = self.fits_image.world_to_image(uniform_names_cat['ra'], uniform_names_cat['dec'], unit='deg')
    uniform_names_cat['x'] = x
    uniform_names_cat['y'] = y
    
    yesno = 'y'
    if colnames_dict['a'] is not None and not self.use_default_names :
        yesno = input("'a', 'b' and 'theta' columns found in catalog. Use them as ellipticity parameters (if not, sources will be shown as circles)? [y] or [n]")
    if colnames_dict['a'] is None or yesno != 'y' :
        size = np.min([self.fits_image.image_data.shape[0], self.fits_image.image_data.shape[1]]) / 1000
        uniform_names_cat['a'] = np.full(len(uniform_names_cat), size)
        uniform_names_cat['b'] = np.full(len(uniform_names_cat), size)
        uniform_names_cat['theta'] = np.full(len(uniform_names_cat), 0.)
    return uniform_names_cat


def initialize_catalog(cat, self) :
    ##### Combine or just open the catalog #####
    path = None
    if isinstance(cat, str) :
        path = cat
        cat = open_cat(cat)
    elif isinstance(cat, list) :
        path = os.path.dirname(cat[0])
        run_match(cat[0], cat[1])
        for i in range(len(cat)-3) :
            run_match('matched_A_B.fits', cat[i+2])
        matched_cat = run_match('matched_A_B.fits', cat[-1])
        cat = Table(matched_cat[1].data)
        os.remove('matched_A_B.fits')
    
    ##### Standardize the column names #####
    uniform_names_cat = make_uniform_names_cat(cat, self)
    return uniform_names_cat, path
    







class catalog :
    def __init__(self, cat, fits_image, color=[0., 1., 1., 0., 0.5], use_default_names=True, units=None, verbose=True) :
        self.fits_image = fits_image
        
        self.xy_axes = None
        self.use_default_names = use_default_names
        self.units = units
        self.verbose = verbose
        
        self.cat, self.path = initialize_catalog(cat, self)
        self.qtItems = [] #np.empty(len(self.cat), dtype=PyQt5.QtWidgets.QGraphicsEllipseItem)
        #self.qtItems = [PyQt5.QtWidgets.QGraphicsEllipseItem() for _ in range(len(self.cat))]
        self.qtItems_column = [] #np.empty(len(self.cat), dtype=pg.TextItem)
        #self.qtItems = np.empty(len(self.cat), dtype=utils.utils_classes.selectable_ellipse.SelectableEllipse)
        self.color = color if color is not None else [0., 1., 1., 0., 0.5]
        self.selection_mask = np.full(len(self.cat), False)
        self.selection_regions = []
        self.Scatter_widget = None
        self.x_axis_cleaned = np.full(len(self.cat), None)
        self.y_axis_cleaned = np.full(len(self.cat), None)
        
        #self.is_plotted = False
        #self.is_plotted_column = False
    
    
    def _vprint(self, *args, **kwargs) :
        if self.verbose :
            print(*args, **kwargs)
    
    def make_mask_naninf(self, xy_axes) :
        x_axis = self.cat[xy_axes[0]]
        y_axis = self.cat[xy_axes[1]]
        
        nan_mask = np.logical_not(np.isnan(x_axis)) & np.logical_not(np.isnan(y_axis))
        inf_mask = (x_axis!=np.inf) & (y_axis!=np.inf)
        #extremes_mask = (x_axis>0) & (x_axis<50)
        self.mask_naninf = nan_mask & inf_mask
        self.x_axis_cleaned = x_axis[self.mask_naninf]
        self.y_axis_cleaned = y_axis[self.mask_naninf]
        
        self.cat = self.cat[self.mask_naninf]
        self.plot()
        #self.qtItems = np.empty(len(self.cat), dtype=PyQt5.QtWidgets.QGraphicsEllipseItem)
        self.selection_mask = np.full(len(self.cat), False)
    
    def plot(self, scale=1., color=None, text_column=None, linewidth=3, marker=None) :
        self.clear()
        self._vprint('Plotting...')
        x = self.cat['x']
        y = self.cat['y']
        semi_major = self.cat['a'] * scale
        semi_minor = self.cat['b'] * scale
        angle = self.cat['theta']
        for i in tqdm(range(len(semi_major))) :
            ellipse = self.plot_one_object(x[i], y[i], semi_major[i], semi_minor[i], angle[i], i, \
                                           color=color, linewidth=linewidth, marker=marker, size=scale*15.)
            self.qtItems.append(ellipse)
        
        if text_column is not None:
            self.plot_column(text_column, color=color)
    
    def plot_column(self, text_column, color=None, n_digit=3):
        self.clear_column()
        self._vprint('Plotting...')
        if text_column not in self.cat.colnames:
            self._vprint(f"Column '{text_column}' not found in catalog")
            return
        
        color = self.color if color is None else color
        color = list(np.array(self.color[:3])*255)
        
        for i in tqdm(range(len(self.cat))) :
            to_plot = self.cat[text_column][i]
            text = to_plot if type(to_plot)==str else f"{to_plot:.{n_digit}g}"
            text_item = pg.TextItem(text, color=color)
    
            # Get ellipse position and size for offset calculation
            x = self.cat['x'][i]
            y = self.fits_image.image_data.shape[0] - self.cat['y'][i]  # Flip y to match PyQtGraph convention
            semi_major = self.cat['a'][i]
            semi_minor = self.cat['b'][i]
    
            # Position text slightly offset from the ellipse
            offset = max(semi_major, semi_minor)
            text_item.setPos(x + offset/2, y - offset/2)
    
            font = PyQt5.QtGui.QFont()
            font.setPointSize(15)
            text_item.setFont(font)

            self.fits_image.qt_image.addItem(text_item)
            self.qtItems_column.append(text_item)
    
    def clear(self) :
        #qtItems_list = self.fits_image.qt_image.getView().allChildItems()
        self._vprint('Clearing galaxies...')
        for qtItem in self.qtItems :
            self.fits_image.qt_image.removeItem(qtItem)
            del qtItem
        self.qtItems.clear()
        
        self.clear_column()
    
    def clear_column(self) :
        self._vprint('Clearing column labels...')
        for text_item in self.qtItems_column:
            self.fits_image.qt_image.removeItem(text_item)
            del text_item
        self.qtItems_column.clear()
    
    def clear_selection(self) :
        self.selection_mask[np.full(len(self.cat), True)] = False
        self.selection_regions.clear()
        self.clear()
        self.plot()
    
    def plot_one_object(self, x, y, semi_major, semi_minor, angle, idx, color=None, linewidth=3, marker=None, size=15) :
        if color is None :
            color = self.color
        #make the flip to accomosate pyqtgraph's plotting conventions
        y = self.fits_image.image_data.shape[0] - y
        angle = -angle
        #####################################################################
        if marker==None or marker=='ellipse' :
            if self.Scatter_widget==None :
                scatter_pos = None
            else :
                scatter_pos = (self.x_axis_cleaned[idx], self.y_axis_cleaned[idx])
            to_plot = SelectableEllipse(x-semi_major/2, y-semi_minor/2, semi_major, semi_minor, idx, self.selection_mask, 
                                        color, linewidth=linewidth, scatter_pos=scatter_pos, 
                                        Scatter_widget=self.Scatter_widget)
            to_plot.setTransformOriginPoint( PyQt5.QtCore.QPointF(x, y) )
            to_plot.setRotation(angle)
        else :
            to_plot = pg.ScatterPlotItem(size=size, symbol=marker)
            color = make_full_color(color)
            if color[3]==0 : #filled_markers==False
                to_plot.setPen( pg.mkPen(color[:3], width=2) ) #no outline
                to_plot.setBrush( pg.mkBrush([0,0,0,0]) )
            else :
                to_plot.setPen( pg.mkPen([0,0,0,0], width=0.1) ) #outline makes cross thicker
                to_plot.setBrush( pg.mkBrush(color[:-2]) )
                
            to_plot.setData([x], [y])
        
        self.fits_image.qt_image.addItem(to_plot)
        return to_plot
    
    def make_selection_panel(self, xy_axes=None) :
        if xy_axes is None :
            xy_axes = []
            xy_axes.append( input('select x axis among: ' + str(self.cat.colnames)) )
            xy_axes.append( input('select y axis among: ' + str(self.cat.colnames)) )
        else :
            x_axis = self.cat[xy_axes[0]]
            y_axis = self.cat[xy_axes[0]]
        self.xy_axes = xy_axes
        
        self.make_mask_naninf(xy_axes)
        
        self.Scatter_widget = SelectableScatter(self)
        #self.Scatter_widget.setTitle('Red sequence')
        
        self.Scatter_widget.setAspectLocked(lock=True, ratio=1)
        self.Scatter_widget.autoRange()
        #self.Scatter_widget.setSizePolicy(pg.QtWidgets.QSizePolicy.Fixed, pg.QtWidgets.QSizePolicy.Expanding)
        self.fits_image.qt_layout.addWidget(self.Scatter_widget)
        
        if xy_axes is not None :
            self.Scatter_widget.setLabel('bottom', xy_axes[0])
            self.Scatter_widget.setLabel('left', xy_axes[1])
        
        #self.selection_mask = np.full(len(self.mag_F444W_cleaned), False)
        self.plot()
    
    def remove_selection_panel(self) :
        self.Scatter_widget.hide()
        del self.Scatter_widget
        self.Scatter_widget = None
    
    def make_image_ROI(self) :
        """
        Creates an ellipse ROI to add elliptic sources to the catalog by drawing them by hand.
        """
        center_y = self.fits_image.image_data.shape[0]/2
        center_x = self.fits_image.image_data.shape[1]/2
        self.image_ROI = ellipse_maker_ROI([center_x-200, center_y-100], [400, 200], self.fits_image.qt_image, self.fits_image.window, self.cat)
        make_handles(self.image_ROI)
        self.fits_image.qt_image.addItem(self.image_ROI)
        
    def make_selection_ROI(self) :
        """
        Creates a rectangle ROI to select all sources inside of it.
        """
        self.fits_image.image_widget.cat = self
        self.select_sources = SelectSources(self)
        
    def save_selection_mask(self, path=None) :
        self.selection_mask_path = self.make_path(path, self.path, 'selection_mask.npy')
        np.save(self.selection_mask_path, self.selection_mask)
        self._vprint("Selection mask saved at " + self.selection_mask_path)
        
    def load_selection_mask(self, path=None) :
        self.selection_mask_path = self.make_path(path, self.path, 'selection_mask.npy')
        self.selection_mask = np.load(self.selection_mask_path)
        
    def save_selection_regions(self, path=None) :
        self.selection_regions_path = self.make_path(path, self.fits_image.image_path, 'selection_regions.npy')
        np.save(self.selection_regions_path, self.selection_regions)
        
    def load_selection_regions(self, path=None, name='selection_regions.npy') :
        self.selection_regions_path = self.make_path(path, self.fits_image.image_path, name)
        self.selection_regions = np.load(self.selection_regions_path).tolist()
        
        size_y = self.fits_image.qt_image.image.shape[0]
        for rect_params in self.selection_regions :
            indiv_mask = InRectangle(self.cat['x'], size_y - self.cat['y'], rect_params)
            self.selection_mask[indiv_mask] = True
        
        fig, ax = plt.subplots()
        ax.axis('equal')
        size = max(self.fits_image.qt_image.image.shape[0], self.fits_image.qt_image.image.shape[1])
        #ax.invert_yaxis()
        ax.set_ylim([size+2000, -2000])
        ax.set_xlim([-4000, size+4000])
        
        for i, rect_params in enumerate(self.selection_regions) :
            x0, y0, a, b, angle = rect_params
            #x1, y1 = x0, y0
            #x2, y2 = x0 + a*np.cos(angle), y0 + a*np.sin(angle)
            #x3, y3 = x0 + a*np.cos(angle) - b*np.sin(angle),  y0 + a*np.sin(angle) + b*np.cos(angle)
            #x4, y4 = x0 - b*np.sin(angle), y0 + b*np.cos(angle)
            #ax.plot([x1, x2, x3, x4, x1], size_y-np.array([y1, y2, y3, y4, y1]), c='b')
            #ax.axis('equal')
            ax.add_patch( Rectangle((x0, y0), a, b, angle=angle*180/np.pi, alpha=0.4 ))
            ax.text(x0, y0, str(i))
            fig.show()
            plt.pause(0.05)
            
    def make_path(self, path, ref_path, name) :
        if path is None :#and self.ref_path is not None :
            #to_return = os.path.join(os.path.dirname(ref_path), name)
            to_return = os.path.join(os.path.dirname(ref_path), os.path.basename(ref_path).split('.')[0] + '_' + name)
        elif os.path.isdir(path) :
            to_return = os.path.join(path, name)
        elif os.path.isdir(os.path.dirname(path)) :
            to_return = path
        return to_return
        
        
    def plot_one_galaxy_mpl(self, x, y, a, b, theta, color=[1,1,1], text=None, ax=None, linewidth=1., text_color='white', text_alpha=0.5) :
        edgecolor = list(color).copy()
        edgecolor.append(1)
        facecolor = edgecolor.copy()
        facecolor[-1] = 0
        ellipse = Ellipse( (x, y), a, b, angle=theta, facecolor=facecolor, edgecolor=edgecolor, lw=linewidth )
        if ax is None :
            self.fits_image.mpl_ax.add_artist(ellipse)
            if text is not None :
                #self.fits_image.mpl_ax.text(x-1.5*b*np.abs(np.sin(theta)), y-1.5*b*np.abs(np.cos(theta)), text, color=edgecolor[:3], \
                #                 ha='right', va='top')
                offset = 0.85
                theta_modulo = theta%180 * np.pi/180
                if theta_modulo<np.pi/2 :
                    x_text, y_text = x+offset*b*np.abs(np.sin(theta_modulo)), y-offset*b*np.abs(np.cos(theta_modulo))
                    horizontalalignment, verticalalignment = 'left', 'top'
                else :
                    x_text, y_text = x-offset*b*np.abs(np.sin(theta_modulo)), y-offset*b*np.abs(np.cos(theta_modulo))
                    horizontalalignment, verticalalignment = 'right', 'top'
                self.fits_image.mpl_ax.text( x_text, y_text, text, c=text_color, alpha=1, fontsize=15, \
                                              ha=horizontalalignment, va=verticalalignment, \
                                              bbox=dict(facecolor=edgecolor[:3], alpha=text_alpha, edgecolor='none') )
        else :
            ax.add_artist(ellipse)
            if text is not None :
                offset = 0.85
                theta_modulo = theta%180 * np.pi/180
                if theta_modulo<np.pi/2 :
                    x_text, y_text = x+offset*b*np.abs(np.sin(theta_modulo)), y-offset*b*np.abs(np.cos(theta_modulo))
                    horizontalalignment, verticalalignment = 'left', 'top'
                else :
                    x_text, y_text = x-offset*b*np.abs(np.sin(theta_modulo)), y-offset*b*np.abs(np.cos(theta_modulo))
                    horizontalalignment, verticalalignment = 'right', 'top'
                #ax.text(x_text, y_text, text, color=edgecolor[:3], \
                #        ha=horizontalalignment, va=verticalalignment)
                ax.text( x_text, y_text, text, c=text_color, alpha=1, fontsize=15, \
                         ha=horizontalalignment, va=verticalalignment, \
                         bbox=dict(facecolor=edgecolor[:3], alpha=text_alpha, edgecolor='none') )
    
    def export_to_mult_file(self, file_path=None) :
        if file_path is None :
            file_path = os.path.join(os.path.dirname(self.fits_image.image_path), 'mult.lenstool')
        
        sub_cat = self.cat[self.selection_mask] if True in self.selection_mask else self.cat
        
        header = "#REFERENCE 0\n## id   RA      Dec        a         b         theta     z         mag\n"
        
        
        theta_colname, z_colname, mag_colname = None, None, None
        
        # Check for theta column
        if 'THETA_WORLD' in sub_cat.colnames :
            theta_colname = 'THETA_WORLD'
        else :
            for name in sub_cat.colnames :
                if 'theta' in name.lower() :
                    theta_colname = name
                    print(f"Using {theta_colname} column in exported catalog, make sure {theta_colname} is relative to world coordinates and not image coordinates!")
                    break
        
        # Check for redshift colum
        for name in ['z', 'z_spec', 'zspec', 'z_phot', 'zphot', 'zb', 'redshift'] :
            colnames_lower = [name.lower() for name in sub_cat.colnames]
            if name in colnames_lower :
                z_colname = sub_cat.colnames[colnames_lower.index(name)]
                break
        
        # Check for magnitude column
        for name in sub_cat.colnames :
            if 'mag' in name.lower() :
                mag_colname = name
                break
        
        with open(file_path, 'w') as f :
            f.write(header)
            for index, row in enumerate(sub_cat) :
                theta = row[theta_colname] if theta_colname is not None else 0.0
                z = row[z_colname] if z_colname is not None else 0.0
                mag = row[mag_colname] if mag_colname is not None else 0.0
                line = (f"{row['id']:<3}  {row['ra']:10.6f}  {row['dec']:10.6f}  "
                        f"{row['a']:8.6f}  {row['b']:8.6f}  {theta:8.6f}  "
                        f"{z:8.6f}  {mag:8.6f}\n")
                f.write(line)
        self._vprint('Selected sources exported at ' + file_path)
                
    
    def export_to_potfile(self, file_path=None, units='pixel') :
        cat = self.cat[self.selection_mask] if True in self.selection_mask else self.cat
        cat = cat.copy()
        
        if units=='pixel' :
            cat['a'] *= self.fits_image.pix_deg_scale*3600
            cat['b'] *= self.fits_image.pix_deg_scale*3600
            self._vprint("Converting pixel units to arcsec")
        elif units=='deg' :
            cat['a'] *= 3600
            cat['b'] *= 3600
            self._vprint("Converting deg units to arcsec")
        else :
            if units!='arcsec' :
                self._vprint("Units not recognized, exporting as is")
        
        mag_col = self.xy_axes[0] if self.xy_axes is not None else input('select magnitude column to be used in Lenstool potfile (press return directly to just populate with zeros): ' + str(self.cat.colnames))
        if mag_col in self.cat.colnames :
            self._vprint(f"Using '{mag_col}' as mag column")
            sort_array = np.argsort(cat[mag_col])
            sorted_cat = cat[sort_array]
        else :
            mag_col = None
            self._vprint("No valid column selected. Creating potfile without magnitude information.")
            sorted_cat = cat.copy()
        
        lines = []
        lines.append('#REFERENCE 0\n')
        lines.append('## id   RA   Dec        a        b        theta     mag       lum\n')
        
        for i, galaxy in enumerate(sorted_cat) :
            mag = galaxy[mag_col] if mag_col is not None else 0.0
            lines.append( '%d %f %f %f %f %f %f 0.\n' % (i+1, galaxy['ra'], galaxy['dec'], galaxy['a'], galaxy['b'], galaxy['theta']-self.fits_image.orientation, mag) )
        
        if file_path is None :
            if self.path is not None :
                file_path = os.path.join( os.path.dirname(self.path), 'exported_potfile.lenstool')
            else :
                file_path = os.path.join( os.path.dirname(self.fits_image.image_path), 'exported_potfile.lenstool')
        self._vprint('Exporting selected sources to ' + file_path)
        with open(file_path, 'w') as file:
            file.writelines(lines)
    
    #def transfer_col(self, col_to_transfer) :
    #    if self.fits_image.imported_cat is not None :
    #        if col_to_transfer in self.fits_image.imported_cat.cat.colnames :
    #            temp_cat = match_cat2([self.cat, self.fits_image.imported_cat.cat], keep_all_col=True, fill_in_value=-1)
    #            if col_to_transfer in self.cat.colnames :
    #                col_to_transfer = col_to_transfer + '_CAT2'
    #            self.cat[col_to_transfer] = temp_cat[col_to_transfer]
    #            print('###############\nColumn ' + col_to_transfer + ' added.\n###############')
    #        else :
    #            print(col_to_transfer + ' not found in imported_cat')
    #    else :
    #        print('No imported_cat')
    
    def transfer_col(self, col_to_transfer, which_cat="imported_cat", index=None):
        if index is not None :
            source_cat = self.fits_image.imported_cat_list[index]
        else :
            source_cat = getattr(self.fits_image, which_cat, None)
    
        if source_cat is not None:
            if col_to_transfer in source_cat.cat.colnames:
                temp_cat, match_idx = match_cat2([self.cat, source_cat.cat], keep_all_col=True, return_match_idx=True) #, fill_in_value=-1.0
                if col_to_transfer in self.cat.colnames:
                    col_to_transfer = col_to_transfer + '_CAT2'
                self.cat[col_to_transfer] = temp_cat[col_to_transfer]
                self._vprint(f'###############\nColumn {col_to_transfer} added.\n###############')
            else:
                self._vprint(f'{col_to_transfer} not found in {which_cat}')
        else:
            self._vprint(f'No {which_cat}')
        return match_idx
    
    def export_to_LaTeX(self, columns=[], file_path=None, use_best_fit_value=False) :
        if len(columns)==0 :
            columns = list(self.cat.colnames)
        
        # Detect triplets/quadruplets involving _16/_50/_84_percentile columns and fold
        # them into a single column with asymmetric uncertainties as super/subscript.
        #
        # use_best_fit_value=True  → value source is the base column (name);
        #                            _50_percentile columns are suppressed entirely.
        # use_best_fit_value=False → value source is name_50_percentile when available,
        #                            falling back to the base column otherwise.
        columns_set = set(columns)
        percentile_map = {}   # base_col -> (col_16, col_value, col_84)
        percentile_cols = set()
        for col in columns :
            if col.endswith('_16_percentile') :
                base   = col[:-len('_16_percentile')]
                col_84 = base + '_84_percentile'
                col_50 = base + '_50_percentile'
                if col_84 not in columns_set :
                    continue
                # Require the base column when using best-fit values; otherwise it is
                # optional because _50_percentile can stand in as the display column.
                if use_best_fit_value and base not in columns_set :
                    continue
                if not use_best_fit_value and base not in columns_set and col_50 not in columns_set :
                    continue

                col_value = base if (use_best_fit_value or col_50 not in columns_set) else col_50
                percentile_map[base] = (col, col_value, col_84)
                percentile_cols.add(col)
                percentile_cols.add(col_84)
                if col_50 in columns_set :
                    percentile_cols.add(col_50)

        # When using the median as central value, the base column is not necessarily
        # in the catalog; add a sentinel so display_columns still contains it.
        if not use_best_fit_value :
            for base, (col_16, col_value, col_84) in percentile_map.items() :
                if base not in columns_set :
                    columns_set.add(base)
                    # Insert the synthetic base name right before its _16_percentile column
                    idx = columns.index(col_16)
                    columns.insert(idx, base)

        # Merge z_in / z_opt into a single 'z' column.
        # z_in takes priority; z_opt is used (marked with *) when z_in is absent or invalid.
        merge_z = 'z_in' in columns_set and 'z_opt' in columns_set
        suppressed_cols = percentile_cols.copy()
        if merge_z :
            suppressed_cols.add('z_opt')

        # Only display the base columns; percentile partners are merged into them.
        display_columns = [col for col in columns if col not in suppressed_cols]

        def _colhead(col) :
            return 'z' if (merge_z and col == 'z_in') else col

        _UNITS_MAP = {
            'ra'    : 'deg',
            'dec'   : 'deg',
            'theta' : 'deg',
            'x'     : 'pix',
            'y'     : 'pix',
            'a'     : 'pix',
            'b'     : 'pix',
        }
        def _units(col) :
            name = _colhead(col)
            if name in _UNITS_MAP :
                return _UNITS_MAP[name]
            if 'mag' in name.lower() :
                return 'mag'
            return ''

        # Build both header lines; only emit the units line when at least one column
        # has a known unit (otherwise the second line would be all empty colheads).
        col_names_row = " &\n".join(f"\\colhead{{{_colhead(col)}}}" for col in display_columns)
        units_entries = [f"({_units(col)})" if _units(col) else "" for col in display_columns]
        has_units     = any(u for u in units_entries)
        units_row     = " &\n".join(f"\\colhead{{{u}}}" for u in units_entries)

        table_str = ("\\begin{deluxetable*}{l" + " c"*len(display_columns) + "}\n"
                     "\\tablecaption{\\label{tab:}}\n"
                     "\\tablewidth{\\columnwidth}\n"
                     "\\tablehead{\n")
        if has_units :
            table_str += col_names_row + " \\\\\n"
            table_str += units_row + "\n"
        else :
            table_str += col_names_row + "\n"
        table_str += ("\\vspace{-0.07in}\\\\}\n"
                      "\\startdata\n")
        for src in self.cat :
            for col in display_columns :
                if merge_z and col == 'z_in' :
                    z_in_val = src['z_in']
                    if not (np.isnan(float(z_in_val)) or float(z_in_val) <= 0) :
                        to_parse = f"{z_in_val:.3g}"
                    else :
                        if 'z_opt' in percentile_map :
                            col_16, col_value, col_84 = percentile_map['z_opt']
                            val    = src[col_value]
                            upper  = src[col_84] - val
                            lower  = val - src[col_16]
                            to_parse = f"${val:.3g}^{{+{upper:.3g}}}_{{-{lower:.3g}}}$"
                        else :
                            z_opt_val = src['z_opt']
                            to_parse = f"{z_opt_val:.3g}*"
                elif col in percentile_map :
                    col_16, col_value, col_84 = percentile_map[col]
                    val    = src[col_value]
                    upper  = src[col_84] - val
                    lower  = val - src[col_16]
                    to_parse = f"${val:.3g}^{{+{upper:.3g}}}_{{-{lower:.3g}}}$"
                elif type(src[col]) in [np.str_, str] :
                    to_parse = src[col]
                elif col in ['ra', 'dec'] :
                    to_parse = f"{src[col]:#.8g}"
                else :
                    to_parse = f"{src[col]:.3g}"
                table_str += to_parse + " &\n"
            table_str = table_str[:-3] + "\\\\\n"
        table_str += ("\\enddata\n"
                      "\\tablecomments{}\n"
                      "\\end{deluxetable*}")
        self._vprint(table_str)
        
        if file_path is None and self.path is not None :
            file_path = os.path.join(os.path.dirname(self.path), 'catalog_latex_table.txt')
        if file_path is not None :
            with open(file_path, 'w') as f :
                f.write(table_str)
            self._vprint('LaTeX table saved at ' + file_path)
        return table_str
    
    #def export_thumbnails(self, mask=None, group_images=True) :
    #    if mask is None :
    #        mask = self.selection_mask
    #    print('Function is being built.')
    #    ### TO DO ###







