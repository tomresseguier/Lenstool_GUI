import numpy as np
import os
import glob
import pyqtgraph as pg
from astropy.table import Table
from collections import defaultdict
from PyQt5.QtCore import QRectF

from .utils_multiple_images import import_multiple_images
from ...utils.utils_astro.utils_general import relative_to_world
from ...utils.utils_general.utils_general import find_close_coord



def get_lenstool_file_path(model_dir, name) :
    path_list = glob.glob( os.path.join(model_dir, f"*{name}*.lenstool") )
    if len(path_list)==0 :
        return None
    else :
        if len(path_list)==1 :
            print(f"{name} file found: '" + path_list[0] + "'")
        elif len(path_list)>1 :
            print(f"Several {name} files found: ")
            print(path_list)
            print("Using first item in list: '" + path_list[0] + "'")
        return path_list[0]


def import_sources(self, predicted_sources_path, image, AttrName='source', units='pixel', filled_markers=False) :
    with open(predicted_sources_path) as file :
        source_lines = file.readlines()[1:]
    sources = Table(names=['id','ra','dec','a','b','theta','z','mag'], dtype=['str', *['float',]*7])
    for line in source_lines :
        sources.add_row(line.split())
    sources['ra'], sources['dec'] = self.relative_to_world(sources['ra'], sources['dec'])
    self.source = image.make_catalog(sources, color=[1.,0.5,0.], units='arcsec', verbose=self.verbose)


def export_thumbnails(self, group_images=True, square_thumbnails=True, square_size=150, margin=50, distance=200, export_dir=None, boost=True, make_broad_view=True, broad_view_params=None) :
    export_dir = os.path.join(os.path.dirname(self.image.image_path), 'mult_thumbnails') if export_dir is None else os.path.abspath(os.path.join(export_dir, 'mult_thumbnails'))
    if os.path.isdir( os.path.dirname( os.path.dirname(export_dir) ) ) and not os.path.isdir( os.path.dirname(export_dir) ) :
        os.mkdir(os.path.dirname(export_dir))
    if not os.path.isdir(export_dir) :
        os.mkdir(export_dir)
    
    if not self.image.boosted and boost :
        self.image.boost()
    
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
        
        #if square_thumbnails :
        x_side_size = x_max - x_min
        y_side_size = y_max - y_min
        demi_taille_unique = round( max(x_side_size, y_side_size)/2 )
        if x_side_size!=y_side_size :
            x_pix = round( (x_max + x_min)/2 )
            y_pix = round( (y_max + y_min)/2 )
            x_min = x_pix - demi_taille_unique
            x_max = x_pix + demi_taille_unique
            y_min = y_pix - demi_taille_unique
            y_max = y_pix + demi_taille_unique
        
        zoom_rect = QRectF(x_min, self.image.image_data.shape[0] - y_max, demi_taille_unique*2, demi_taille_unique*2)
        self.image.ImageView.getView().setRange(zoom_rect)        
        
        
        thumbnail_path = os.path.join( export_dir, 'mult_' + group[0] + '.png' )
        print('Creating ' + thumbnail_path)
        exporter = pg.exporters.ImageExporter(self.image.ImageView.view)
        exporter.export(thumbnail_path)
        print('Done')
        
    
    ##### Adding broad view #####
    if make_broad_view :
        # broad_view_params = [ [x_min, x_max], [y_min, y_max] ]
        if broad_view_params is not None :
            x = broad_view_params[0][0]
            y = self.image.image_data.shape[0] - broad_view_params[1][1]
            x_width = broad_view_params[0][1]-broad_view_params[0][0]
            y_width = broad_view_params[1][1]-broad_view_params[1][0]
            zoom_rect = QRectF(x, y, x_width, y_width)
        else :
            zoom_rect = QRectF(0, 0, self.image.image_data.shape[1], self.image.image_data.shape[0])
        self.image.ImageView.getView().setRange(zoom_rect)        
        
        broadview_filename = 'broadview'
        for name in self.which :
            broadview_filename += '_' + name
        broadview_path = os.path.join( export_dir, broadview_filename + '.png' )
        print('Creating ' + broadview_path)
        exporter = pg.exporters.ImageExporter(self.image.ImageView.view)
        exporter.export(broadview_path)
        print('Done')
    ##############################


def add_optimized_redshifts(mult, param_best) :
    if not 'z_opt' in mult.colnames :
        mult.add_column(np.full(len(mult), np.nan), name='z_opt')
    if param_best is not None :
        if 'image' in param_best :
            if 'z_m_limit' in param_best['image'] :
                for l in param_best['image']['z_m_limit'] :
                    if type(l[1]) is str :
                        names = []
                        i = 1
                        while type(l[i]) is str :
                            names.append(l[i])
                            i+=1
                    else :
                        names = [str(l[1])]
                        i=2
                    
                    z = l[i+1]
                    for name in names :
                        if name in mult['id'] :
                            fam = mult[mult['id']==name][0]['family']
                            mult['z_opt'][mult['family']==fam] = z


class curves :
    def __init__(self, curves_dir, LensModel, image, which_critcaus='critical', join=False, size=2) :
        self.dir = curves_dir
        self.paths = glob.glob(os.path.join(curves_dir, "*.dat"))
        self.LensModel = LensModel
        self.image = image
        self.size = size
        
        self.qtItems = {}
        for name in LensModel.broad_families :
            self.qtItems[name] = None
        
        self.coords = {}
        for name in LensModel.broad_families :
            curve_mask = np.array([name in os.path.basename(path) for path in self.paths])
            if True in curve_mask :
                lines = []
                for curve_path in np.array(self.paths)[curve_mask] :
                    file = open(curve_path, 'r')
                    all_lines = file.readlines()
                    lines += all_lines[1:]
                
                ra_ref = float( all_lines[0].split()[-2] )
                dec_ref = float( all_lines[0].split()[-1] )
                
                if which_critcaus=='critical' :
                    delta_ra = np.array( [ float( lines[i].split()[1] ) for i in range(len(lines)) ] )
                    delta_dec = np.array( [ float( lines[i].split()[2] ) for i in range(len(lines)) ] )
                    ra, dec = relative_to_world(delta_ra, delta_dec, (ra_ref, dec_ref))
                    x, y = image.world_to_image(ra, dec)
                    
                    y = image.image_data.shape[0] - y
                    
                    shorten_indices = np.linspace(0, len(x) - 1, 10000, dtype=int)
                    x = x[shorten_indices]
                    y = y[shorten_indices]
                    
                    #if join :
                    #    x, y = rearrange_points(x, y)
                         
                if which_critcaus=='caustic' :
                    delta_ra = np.array( [ float( lines[i].split()[3] ) for i in range(len(lines)) ] )
                    delta_dec = np.array( [ float( lines[i].split()[4] ) for i in range(len(lines)) ] )
                    ra, dec = relative_to_world(delta_ra, delta_dec, (ra_ref, dec_ref))
                    x, y = image.world_to_image(ra, dec)
                    
                    y = image.image_data.shape[0] - y
                    
                    shorten_indices = np.linspace(0, len(x) - 1, 10000, dtype=int)
                    x = x[shorten_indices]
                    y = y[shorten_indices]
                    
                    #if join :
                    #    x, y = rearrange_points(x, y)
                
                self.coords[name] = (x, y)
        
    def plot(self) :
        self.clear()
        for name in self.LensModel.which :
            color = np.round(self.LensModel.mult_colors(saturation=self.LensModel.saturation)[name]*255).astype(int)
            color[3] = 255
            
            x, y = self.coords[name]
            
            scatter = pg.ScatterPlotItem(x, y, pen=None, brush=pg.mkBrush(color), size=self.size)
            self.qtItems[name] = scatter
            self.image.ImageView.addItem(scatter)
            
    def clear(self) :
        for name, qtItem in self.qtItems.items() :
            if qtItem is not None :
                self.image.ImageView.removeItem(qtItem)
                self.qtItems[name] = None
                
    








def find_families(image_ids):
    family_ids = image_ids.copy()
    confidence = np.full(len(family_ids), 2)
    for i, name in enumerate(family_ids) :
        if name.startswith('cc') :
            family_ids[i] = name[2:]
            confidence[i] = 0
        elif name.startswith('c') :
            family_ids[i] = name[1:]
            confidence[i] = 1
    
    families = find_families_part2(family_ids)
    
    combined_families = families.copy()
    for i, family in enumerate(families) :
        prefix1 = family.split('.')[0]
        for fam in families :
            prefix2 = fam.split('.')[0]
            if prefix1==prefix2 and family!=fam :
                combined_families[i] = prefix1 + '.'
    #combined_families = np.unique(combined_families)
    
    broad_families = combined_families.copy()
    letter_id = []
    for i, family in enumerate(combined_families) :
        if family[0].isalpha() :
            letter_id.append(i)
            for fam in combined_families :
                if fam[0]==family[0] :
                    broad_families[i] = family[0]
    
    
    
    families_int = families.copy()
    for i in letter_id :
        families_int[i] = str( ord( families[i][0].lower() )-96 )
    
    families_sorted, indices = np.unique(families, return_index=True)
    families_sorted_int = np.array(families_int)[indices]
    
    families_sorted_int = [int(family.split('.')[0]) for family in families_sorted_int]
    families_sorted = families_sorted[np.argsort(families_sorted_int)]
    
    
    
    broad_families_int = broad_families.copy()
    for i in letter_id :
        broad_families_int[i] = str( ord( broad_families[i][0].lower() )-96 )
    
    broad_families_sorted, indices = np.unique(broad_families, return_index=True)
    broad_families_sorted_int = np.array(broad_families_int)[indices]
    
    broad_families_sorted_int = [int(family.split('.')[0]) for family in broad_families_sorted_int]
    broad_families_sorted = broad_families_sorted[np.argsort(broad_families_sorted_int)]
    
    return families, broad_families, families_sorted.tolist(), broad_families_sorted.tolist(), confidence


def find_families_part2(image_ids) :
    # Step 1: Initial guess by chopping last character
    id_to_family = {img_id: img_id[:-1] for img_id in image_ids}
    
    # Step 2: Group by these tentative families
    family_groups = defaultdict(list)
    for img_id, fam in id_to_family.items():
        family_groups[fam].append(img_id)

    # Step 3: Merge singleton families if their name starts with another family name
    updated = True
    while updated:
        updated = False
        singletons = {fam for fam, ids in family_groups.items() if len(ids) == 1}
        for fam in list(singletons):
            for target in family_groups:
                if fam != target and fam.startswith(target):
                    family_groups[target].extend(family_groups[fam])
                    del family_groups[fam]
                    updated = True
                    break
            if updated:
                break

    # Step 4: Merge families with 'alt' in original IDs if the ID starts with another family name
    for fam in list(family_groups):
        for img_id in family_groups[fam]:
            if 'alt' in img_id:
                for target in family_groups:
                    if fam != target and img_id.startswith(target):
                        family_groups[target].extend(family_groups[fam])
                        del family_groups[fam]
                        break
                break  # Only need to check one 'alt' image to trigger a merge

    # Step 5: Build final output mapping
    final_map = {}
    for fam, ids in family_groups.items():
        for img_id in ids:
            final_map[img_id] = fam

    return [final_map[img_id] for img_id in image_ids]




def import_lenstool_files(self) :
    if self.arclets is None :
        arclets_path_list = glob.glob( os.path.join(self.model_dir, "*arclet*.lenstool") )
        if len(arclets_path_list)==1 :
            arclets_path = arclets_path_list[0]
            print(f"{os.path.basename(arclets_path)} found and used as arclets.")
            import_multiple_images(self, arclets_path, self.image, AttrName='arclets', units='pixel', filled_markers=False)
                
    if self.images is None :
        predicted_images_path = os.path.join(self.model_dir, 'image.dat')
        if os.path.isfile(predicted_images_path) :
            import_multiple_images(self, predicted_images_path, self.image, AttrName='images', units='pixel', filled_markers=False)
            import_multiple_images(self, predicted_images_path, self.image, AttrName='images_filtered', units='pixel', filled_markers=False)
            self.filter_image()
    
    if self.curves is None :
        curves_dir = os.path.join(self.model_dir, 'curves')
        if os.path.isdir(curves_dir) :
            self.curves = curves(curves_dir, self, self.image, which_critcaus='critical', join=False, size=2)
    
    if self.source is None and self.reference is not None :
        predicted_sources_path = os.path.join(self.model_dir, 'source.dat')
        if os.path.isfile(predicted_sources_path) :
            import_sources(self, predicted_sources_path, self.image, AttrName='source', units='pixel', filled_markers=False)

