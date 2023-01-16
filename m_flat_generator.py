import os
from astropy.io import fits
import cssproc

def create_master_flat(folder):
    flat_files = [f for f in os.listdir(folder) if f.endswith('.fits') and 'flat' in f]

    # Tworzenie słownika z filtrami i odpowiadającymi im listami plików
    filters = {}
    for flat_file in flat_files:
        hdu = fits.open(os.path.join(folder, flat_file))
        filter_name = hdu[0].header['FILTER']
        if filter_name in filters:
            filters[filter_name].append(flat_file)
        else:
            filters[filter_name] = [flat_file]
        hdu.close()

    for filter_name, flat_list in filters.items():
        master_flat = cssproc.median_combine(flat_list)
        master_flat.writeto(os.path.join(folder, 'master_flat_' + filter_name + '.fits'))

create_master_flat('/path/to/folder')
