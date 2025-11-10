import argparse
from .fits_image import fits_image

def main():
    parser = argparse.ArgumentParser(description='Create a FITS image object from a file path.')
    parser.add_argument('fits_file_path', help='The path to the FITS file.')
    args = parser.parse_args()

    try:
        image = fits_image(args.fits_file_path)
        print(f"Successfully created a fits_image instance from: {args.fits_file_path}")
        # You can now work with the 'image' object.
        # For example, you can access its attributes:
        # print(image.image_data)
    except Exception as e:
        print(f"Error creating fits_image instance: {e}")

if __name__ == '__main__':
    main()
