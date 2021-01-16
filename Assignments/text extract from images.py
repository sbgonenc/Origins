from PIL import Image
import pytesseract

pytesseract.pytesseract.tesseract_cmd = r'<full_path_to_your_tesseract_executable>' #command to execute tesseract

img = Image.open('Zeytinli Rock Fest.jpg') #Like text open, Image open
text = pytesseract.image_to_string(img)    # Identifies str in image