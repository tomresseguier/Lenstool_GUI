import numpy as np
import os
import glob
import types
import matplotlib.pyplot as plt
import PyQt5.QtGui
import pyqtgraph as pg
from astropy.table import Table

from ...utils.utils_plots.plot_utils_general import make_palette, adjust_luminosity, adjust_contrast, plot_scale_bar, plot_image_mpl
from ...utils.utils_astro.utils_general import relative_to_world
from ...utils.utils_astro.cat_manip import match_cat2
from ...utils.utils_general.utils_general import find_close_coord
from ...utils.utils_plots.plt_framework import plt_framework




def make_which_colors(self, filled_markers=False, saturation=None) :
    which = self.broad_families if self.which=='all' else self.which
    
    if filled_markers:
        colors = make_palette(len(which), 1, alpha=0.5, sat_fixed=saturation)
    else:
        colors = make_palette(len(which), 1, alpha=0, sat_fixed=saturation)
    
    which_colors_dict = {}
    for i, name in enumerate(which) :
        which_colors_dict[name] = colors[i]
    #for i, mask in enumerate(to_plot_mask):
    #    for multiple_image in cat[mask]:
    #        which_colors_dict[multiple_image['id']] = colors[i]
    return which_colors_dict


def make_full_color_function(families) :
    n_families = len(families)
    def make_full_color_dict(filled_markers=False, saturation=None) :
        if filled_markers:
            colors = make_palette(n_families, 1, alpha=0.5, sat_fixed=saturation)
        else:
            colors = make_palette(n_families, 1, alpha=0, sat_fixed=saturation)
        full_colors_dict = {}
        for i, family in enumerate(families) :
            full_colors_dict[family] = colors[i]
        return full_colors_dict
    return make_full_color_dict


def import_multiple_images(LensModel, mult_file_path_or_cat, image, units=None, AttrName='mult', filled_markers=False, saturation=0.8) :
    if type(mult_file_path_or_cat)==str :
        multiple_images = Table(names=['id','family','broad_family','ra','dec','a','b','theta','z_in','mag','z_opt', 'z','confidence'], dtype=['str','str','str',*['float',]*10])
        with open(mult_file_path_or_cat, 'r') as mult_file:
            for line in mult_file:
                cleaned_line = line.strip()
                if not cleaned_line.startswith("#") and len(cleaned_line)>0 :
                    split_line = cleaned_line.split()
                    row = [split_line[0], '---', '---'] #split_line[0][:-1]
                    for element in split_line[1:8] :
                        row.append(float(element))
                    row.append(np.nan) # z_opt
                    row.append(np.nan) # z
                    row.append(np.nan) # confidence
                    multiple_images.add_row(row)
        multiple_images['theta'] = multiple_images['theta'] - image.orientation
        no_z_in_mask = multiple_images['z_in']==0.
        multiple_images['z_in'][no_z_in_mask] = np.nan
    else :
        multiple_images = mult_file_path_or_cat.copy()
        
    multiple_images['family'], multiple_images['broad_family'], local_families, local_broad_families, multiple_images['confidence'] = find_families(multiple_images['id'])
    add_optimized_redshifts(multiple_images, LensModel.param_best)
    
    if 'z_in' in multiple_images.colnames and 'z_opt' in multiple_images.colnames :
        for i, m in enumerate(multiple_images) :
            m['z'] = m['z_in'] if not np.isnan(m['z_in']) else m['z_opt']
    
    setattr(LensModel, AttrName, image.make_catalog(multiple_images, units=units, verbose=LensModel.verbose))
    getattr(LensModel, AttrName).workspace = LensModel.workspace
    
    if AttrName=='mult' :
        LensModel.families, indices = np.unique(LensModel.families + local_families, return_index=True)
        LensModel.families = LensModel.families[np.argsort(indices)].tolist()
        
        LensModel.broad_families, indices = np.unique(LensModel.broad_families + local_broad_families, return_index=True)
        LensModel.broad_families = LensModel.broad_families[np.argsort(indices)].tolist()
        
        LensModel.which = LensModel.broad_families.copy()
        LensModel.mult_colors = make_full_color_function(LensModel.broad_families) #make_full_color_function(LensModel.broad_families)
        
    getattr(LensModel, AttrName).saturation = saturation
    
    def make_to_plot_masks() :
        to_plot_masks = {}
        #for i, name in enumerate(LensModel.which) :
        #    other_names_mask = np.full(len(LensModel.which), True)
        #    other_names_mask[i] = False
        #    other_names = np.array(LensModel.which)[other_names_mask]
        #    ambiguous_names = []
        #    for other_name in other_names :
        #        if other_name.startswith(name) :
        #           ambiguous_names.append(other_name)
        #    to_plot_mask = np.full(len(getattr(LensModel, AttrName).cat), False)
        #    for j, im_id in enumerate(getattr(LensModel, AttrName).cat['id']) :
        #        if im_id.startswith(name) and True not in [ im_id.startswith(ambiguous_name) for ambiguous_name in ambiguous_names ] :
        #            to_plot_mask[j] = True
        #    to_plot_masks[name] = to_plot_mask
        for family in LensModel.which :
            to_plot_masks[family] = getattr(LensModel, AttrName).cat['broad_family'] == family
            if len(np.unique(to_plot_masks[family]))==1 and np.unique(to_plot_masks[family])[0]==False :
                to_plot_masks[family] = getattr(LensModel, AttrName).cat['family'] == family
                if len(np.unique(to_plot_masks[family]))==1 and np.unique(to_plot_masks[family])[0]==False :
                    to_plot_masks[family] = getattr(LensModel, AttrName).cat['id'] == family
        return to_plot_masks
    def make_overall_mask() :
        overall_mask = np.full(len(getattr(LensModel, AttrName).cat), False)
        for mask in make_to_plot_masks().values() :
            overall_mask = np.logical_or(overall_mask, mask)
        return overall_mask
    getattr(LensModel, AttrName).masks = make_to_plot_masks
    getattr(LensModel, AttrName).mask = make_overall_mask
        
    def plot_multiple_images(self, scale=1, marker=None, filled_markers=filled_markers, color=None, mpl=False, fontsize=9,
                             make_thumbnails=False, square_size=150, margin=50, distance=200, savefig=False, square_thumbnails=True,
                             boost=[2,1.5,1], linewidth=1.7, text_color='white', text_alpha=0.5, saturation=saturation) :
        self.clear()
        self.saturation = saturation
        
        if color is not None :
            colors_dict = {}
            if type(color[0])==list :
                for i, family in enumerate(LensModel.which) :
                    colors_dict[family] = color[i]
            else :
                for i, family in enumerate(LensModel.which) :
                    colors_dict[family] = color
        else :
            colors_dict = LensModel.mult_colors(filled_markers=filled_markers, saturation=saturation)
        
        cat_contains_ellipse_params = len(np.unique(self.cat['a']))!=1
        count = 0
        for name, mask in self.masks().items() :
            # name might be id, family, or broad_family, so we need to look at the correct column
            for colname in ['id', 'family', 'broad_family'] :
                if name in self.cat[colname] :
                    colname_to_use = colname
            broad_family = self.cat['broad_family'][ np.where(self.cat[colname_to_use]==name)[0][0] ]
            for multiple_image in self.cat[mask] :
                # Remove the *1000
                if not cat_contains_ellipse_params :
                    a, b = scale*40, scale*40
                else :
                    a, b = multiple_image['a']*scale, multiple_image['b']*scale
                color = colors_dict[broad_family].copy()
                if multiple_image['confidence']==1 :
                    color/=2
                elif multiple_image['confidence']==0 :
                    color/=4
                ellipse = self.plot_one_object(multiple_image['x'], multiple_image['y'], a, b, 
                                               multiple_image['theta'], count, color=color, 
                                               linewidth=linewidth, marker=marker, size=scale*15)
                #self.qtItems[count] = ellipse
                self.qtItems.append(ellipse)
                count += 1
                
                if mpl :
                    font = {'size':fontsize, 'family':'DejaVu Sans'}
                    plt.rc('font', **font)
                    self.plot_one_galaxy_mpl(multiple_image['x'], multiple_image['y'], a, b, multiple_image['theta'], color=colors_dict[broad_family][:3],
                                             text=multiple_image['id'], linewidth=linewidth, text_color=text_color, text_alpha=text_alpha)
                    #self.plot_one_galaxy_mpl(multiple_image['x'], multiple_image['y'], a, b, multiple_image['theta'], color=colors_dict[broad_family][:3], text=multiple_image['id'])
        
        
        if make_thumbnails :
            if boost is not None :
                adjusted_image = adjust_contrast(self.image.image_data, boost[0], pivot=boost[1])
                adjusted_image = adjust_luminosity(adjusted_image, boost[2])
            else :
                adjusted_image = self.image.image_data
            
            group_images = True
            if group_images :
                group_list = find_close_coord(self.cat[self.mask()], distance)
            else :
                group_list = [[name] for name in self.cat[self.mask()]['id']]
            
            for group in group_list :
                
                x_array = [self.cat[np.where(self.cat['id']==name)[0][0]]['x'] for name in group]
                y_array = [self.cat[np.where(self.cat['id']==name)[0][0]]['y'] for name in group]
                
                x_pix = (np.max(x_array) + np.min(x_array)) / 2
                y_pix = (np.max(y_array) + np.min(y_array)) / 2
                
                half_side = square_size // 2
                
                x_min = round( max( min( np.min(x_array) - margin, x_pix - half_side ), 0) )
                x_max = round( min( max( np.max(x_array) + margin, x_pix + half_side ), self.image.image_data.shape[1]) )
                y_min = round( max( min( np.min(y_array) - margin, y_pix - half_side ), 0) )
                y_max = round( min( max( np.max(y_array) + margin, y_pix + half_side ), self.image.image_data.shape[0]) )
                
                if square_thumbnails :
                    x_side_size = x_max - x_min
                    y_side_size = y_max - y_min
                    if x_side_size!=y_side_size :
                        demi_taille_unique = round( max(x_side_size, y_side_size)/2 )
                        x_pix = round( (x_max + x_min)/2 )
                        y_pix = round( (y_max + y_min)/2 )
                        x_min = x_pix - demi_taille_unique
                        x_max = x_pix + demi_taille_unique
                        y_min = y_pix - demi_taille_unique
                        y_max = y_pix + demi_taille_unique
                
                plt_framework(image=True, figsize=3, drawscaler=1.2)
                font = {'size':9, 'family':'DejaVu Sans'}
                plt.rc('font', **font)
                
                
                #cropped_image = self.image.image_data[y_min:y_max, x_min:x_max, :]
                fig, ax = plot_image_mpl(adjusted_image, wcs=None, wcs_projection=False, units='pixel',
                                         pos=111, make_axes_labels=False, make_grid=False, crop=[x_min, x_max, y_min, y_max])
                
                for multiple_image_id in group :
                    multiple_image = self.cat[np.where(self.cat['id']==multiple_image_id)[0][0]]
                    color = colors_dict[multiple_image_id[:-1]]
                    if not cat_contains_ellipse_params :
                        a, b = 75, 75
                    else :
                        a, b = multiple_image['a'], multiple_image['b']
                    self.plot_one_galaxy_mpl(multiple_image['x']-x_min, multiple_image['y']-y_min, a, b, multiple_image['theta'],
                                             color=color[:3], text=multiple_image['id'], ax=ax, linewidth=linewidth, text_color=text_color, text_alpha=text_alpha)
                
                ax.axis('off')
                #plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
                
                fig.show()
                    
                plot_scale_bar(ax, deg_per_pix=self.image.pix_deg_scale, unit='arcsec',
                               length=1 , color='white', linewidth=2, text_offset=0.01)
                if savefig :
                    fig.savefig(os.path.join(os.path.dirname(self.image.image_path), 'mult_' + group[0]), bbox_inches='tight', pad_inches=0)
                    
                plt_framework(full_tick_framework=True, ticks='out', image=True, width='full', drawscaler=0.8, tickscaler=0.5, minor_ticks=False)
    
    
    def plot_multiple_images_column(self, text_column, which='all', color=None, bbox=0.2) :
        if text_column not in self.cat.colnames:
            LensModel._vprint(f"Column '{text_column}' not found in catalog")
            return
        if not hasattr(self, 'qtItems_column'):
            self.qtItems_column = []
        for text_item in self.qtItems_column:
            self.image.ImageView.removeItem(text_item)
        self.qtItems_column.clear()
        
        if color is not None :
            colors_dict = {}
            if type(color[0])==list :
                for i, family in enumerate(LensModel.which) :
                    colors_dict[family] = color[i]
            else :
                for i, family in enumerate(LensModel.which) :
                    colors_dict[family] = color
        else :
            colors_dict = LensModel.mult_colors(filled_markers=False, saturation=self.saturation)
        
        colors_dict_background = LensModel.mult_colors(filled_markers=False, saturation=self.saturation)

        for name, mask in self.masks().items() :
            #broad_family = name #self.cat['broad_family'][ np.where(self.cat['family']==name)[0][0] ]
            
            ###
            if name in self.cat['broad_family'][mask] :
                broad_family = name
            else :
                broad_family = self.cat['broad_family'][mask][0]
            ###
            
            for multiple_image in self.cat[mask] :
                
                text = str(multiple_image[text_column])

                background_color = None
                if len(text) > 0 : # To avoid background for empty text
                    if type(bbox) is list or type(bbox) is tuple or type(bbox) is np.ndarray :
                        background_color = np.array(bbox)*255
                    elif type(bbox) is float :
                        background_color = list( np.array(colors_dict_background[broad_family])*255 )[:3] + [bbox*255]
                
                if background_color is not None :
                    text_item = pg.TextItem( text, color=list( np.array(colors_dict[broad_family])*255 )[:3], fill=pg.mkBrush(background_color) )
                else :
                    text_item = pg.TextItem( text, color=list( np.array(colors_dict[broad_family])*255 )[:3] )
                
                x = multiple_image['x']
                y = self.image.image_data.shape[0] - multiple_image['y']  # Flip y to match PyQtGraph convention
                semi_major = multiple_image['a']
                semi_minor = multiple_image['b']
                offset = max(semi_major, semi_minor)
                text_item.setPos(x + offset/2, y - offset/2)
                
                font = PyQt5.QtGui.QFont()
                font.setPointSize(15)
                text_item.setFont(font)
                
                self.image.ImageView.addItem(text_item)
                self.qtItems_column.append(text_item)
    
    
    getattr(LensModel, AttrName).plot = types.MethodType(plot_multiple_images, getattr(LensModel, AttrName))
    getattr(LensModel, AttrName).plot_column = types.MethodType(plot_multiple_images_column, getattr(LensModel, AttrName))
    
    def transfer_ids(self, id_name='id') :
        imported_cat = LensModel.workspace.catalog if LensModel.workspace is not None else None
        if imported_cat is not None :
            if id_name in imported_cat.cat.colnames :
                temp_cat = match_cat2([self.cat, imported_cat.cat], keep_all_col=True, fill_in_value=-1, column_to_transfer=id_name)
                if id_name in self.cat.colnames :
                    id_name = id_name + '_CAT2'
                self.cat[id_name] = temp_cat[id_name]
                LensModel._vprint('###############\nColumn ' + id_name + ' added.\n###############')
            else :
                LensModel._vprint(id_name + ' not found in imported_cat')
        else :
            LensModel._vprint('No imported_cat')
    
    getattr(LensModel, AttrName).transfer_ids = types.MethodType(transfer_ids, getattr(LensModel, AttrName))
        
    
