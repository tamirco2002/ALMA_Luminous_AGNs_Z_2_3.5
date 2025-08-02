import os
import openpyxl
import re

def convert_units(value, unit):
    if unit == "uJy" or unit == "uJy/beam":
        return float(value) * 1e-6  # Convert from microjansky to jansky
    elif unit == "mJy" or unit == "mJy/beam":
        return float(value) * 1e-3  # Convert from millijansky to jansky
    else:
        return float(value)

def convert_axis_units(value, unit):
    if unit == "arcsec":
        return float(value)
    elif unit == "marcsec":
        return float(value) * 1e-3  # Convert from milliarcsec to arcsec
    else:
        return None

def process_text_file(file_path, excel_sheet):
    file_name = os.path.splitext(os.path.basename(file_path))[0]  # Extract filename without extension

    with open(file_path, 'r') as file:
        lines = file.readlines()

    integrated_pattern = re.compile(r"Integrated:\s+(-?\d+(\.\d+)?)\s+\+/-\s+(-?\d+(\.\d+)?)\s+(\w+)")
    peak_pattern = re.compile(r"Peak:\s+(-?\d+(\.\d+)?)\s+\+/-\s+(-?\d+(\.\d+)?)\s+(\w+/(\w+))")
    rms_pattern = re.compile(r"RMS of input image:\s+(-?\d+(\.\d+)?([eE]-?\d+)?)\s+(\w+)")
    position_pattern = re.compile(r"---\sra:\s+([-+]?\d+:\d+:\d+\.\d+)\s+\+\/-\s+(\d+\.\d+)\s+s\s+\((\d+\.\d+)\s+arcsec")
    dec_pattern = re.compile(r"---\sdec:\s+([-+]?\d+\.\d+\.\d+\.\d+)\s+\+\/-\s+(\d+\.\d+)\s+arcsec")
    deconv_pattern = re.compile(r"Image component size \(deconvolved from beam\) ---\n(.+)")
    deconv_pattern = re.compile(
        r"Image component size \(deconvolved from beam\) ---\n\s+--- major axis FWHM:\s+(\d+(\.\d+)?)\s+\+/-\s+(\d+(\.\d+)?)\s+marcsec\n\s+--- minor axis FWHM:\s+(\d+(\.\d+)?)\s+\+/-\s+(\d+(\.\d+)?)\s+marcsec"
    )
    deconv_pattern2 = re.compile(
        r"Image component size \(deconvolved from beam\) ---\n\s+--- major axis FWHM:\s+(\d+(\.\d+)?)\s+\+/-\s+(\d+(\.\d+)?)\s+arcsec\n\s+--- minor axis FWHM:\s+(\d+(\.\d+)?)\s+\+/-\s+(\d+(\.\d+)?)\s+arcsec"
    )
    not_deconv_pattern = re.compile(r"Image component size \(deconvolved from beam\) ---\n(.+)")
    upper_limit_pattern = re.compile(r"It may be as large as (\d+(\.\d+)?) arcsec x (\d+(\.\d+)?) arcsec")

    pixels_pattern = re.compile(r"Number of pixels used in fit: (\d+)")

    integrated_values = integrated_errors = peak_values = peak_errors = rms_value = ra_sec = dec_arcsec =  None
    deconvolved_status = major_axis = minor_axis = None

    for line in lines:
        pixels_match = pixels_pattern.search(line)
        if pixels_match:
            pixels_value = int(pixels_match.group(1))
            break  # Stop if you only want the first match

    for line in lines:
        integrated_match = integrated_pattern.search(line)
        if integrated_match:
            value, _, error, _, unit = integrated_match.groups()
            integrated_values = convert_units(value, unit)
            integrated_errors = convert_units(error, unit)
            break

    for line in lines:
        peak_match = peak_pattern.search(line)
        if peak_match:
            value, _, error, _, unit1, unit2 = peak_match.groups()
            unit = unit1 + "/" + unit2
            peak_values = convert_units(value, unit1)
            peak_errors = convert_units(error, unit1)
            break

    for line in lines:
        rms_match = rms_pattern.search(line)
        if rms_match:
            value, _, _, unit = rms_match.groups()
            rms_value = convert_units(value, unit)
            break

    for idx, line in enumerate(lines):
        position_match = position_pattern.search(line)
        if position_match:
            ra_sec, ra_error, dec_arcsec = position_match.groups()
            ra_sec = f"'{ra_sec}"  # Add a single quote to avoid Excel's misinterpretation
            break

    for line in lines:
        dec_match = dec_pattern.search(line)
        if dec_match:
            dec_parts = dec_match.group(1).split('.')
            dec_colon_format = ':'.join(dec_parts[:-1]) + '.' + dec_parts[-1]
            dec_arcsec, dec_error = dec_match.groups()
            dec_arcsec = f"'{dec_colon_format}"  # Add a single quote to avoid Excel's misinterpretation
            break

    deconv_match = not_deconv_pattern.search(''.join(lines))
    if deconv_match:
        major_axis_error = minor_axis_error = None
        deconv_info = deconv_match.group(0).strip()
        if "Could not deconvolve source from beam" in deconv_info:
            deconvolved_status = "No"
            major_axis = minor_axis = None
        elif "Component is a point source" in deconv_info:
            deconvolved_status = "Upper limit"
            major_axis = minor_axis = 0
            upper_limit_match = upper_limit_pattern.search(''.join(lines))
            if upper_limit_match:
                major_axis = float(upper_limit_match.group(1))
                minor_axis = float(upper_limit_match.group(3))

        else:
            try:
                pat_num = 1
                deconv_match = deconv_pattern.search(''.join(lines))
                deconv_info = deconv_match.group(0).strip()
            except:
                pat_num = 2
                deconv_match = deconv_pattern2.search(''.join(lines))
                deconv_info = deconv_match.group(0).strip()

            deconvolved_status = "Yes"
            if pat_num == 1:
                major_axis_match = re.search(r"major axis FWHM:\s+(\d+(\.\d+)?)\s+\+/-\s+(\d+(\.\d+)?)\s+marcsec",
                                             deconv_info)
                minor_axis_match = re.search(r"minor axis FWHM:\s+(\d+(\.\d+)?)\s+\+/-\s+(\d+(\.\d+)?)\s+marcsec",
                                             deconv_info)

            else:
                # Check for matches in the respective patterns and capture values
                major_axis_match = re.search(r"major axis FWHM:\s+(\d+(\.\d+)?)\s+\+/-\s+(\d+(\.\d+)?)\s+arcsec",
                                             deconv_info)
                minor_axis_match = re.search(r"minor axis FWHM:\s+(\d+(\.\d+)?)\s+\+/-\s+(\d+(\.\d+)?)\s+arcsec",
                                             deconv_info)
            major_axis = minor_axis = None

            # Extract major axis value and error if present
            if major_axis_match:
                major_axis, _, major_axis_error, _ = major_axis_match.groups()
                if pat_num ==1:
                    major_axis = float(major_axis) * 1e-3  # Convert from milliarcsec to arcsec
                    major_axis_error = float(major_axis_error) * 1e-3  # Convert from milliarcsec to arcsec

            # Extract minor axis value and error if present
            if minor_axis_match:
                minor_axis, _, minor_axis_error, _ = minor_axis_match.groups()
                if pat_num ==1:
                    minor_axis = float(minor_axis) * 1e-3  # Convert from milliarcsec to arcsec
                    minor_axis_error = float(minor_axis_error) * 1e-3  # Convert from milliarcsec to arcsec

    excel_sheet.append([file_name, integrated_values, integrated_errors, peak_values, peak_errors,rms_value,
                        ra_sec,dec_arcsec,deconvolved_status,major_axis,major_axis_error,minor_axis,minor_axis_error,pixels_value])


def process_folder(folder_path, output_file):
    wb = openpyxl.Workbook()
    sheet = wb.active
    sheet.title = "Data"

    sheet.append(["File Name", "Integrated Value (Jy)", "Integrated Error (Jy)", "Peak Value (Jy/beam)",
                  "Peak Error (Jy/beam)", "RMS Value (Jy)","RA (sec)", "DEC (arcsec)",
                  "Deconvolved Status","Major Axis (arcsec)","Major Axis Error (arcsec)","Minor Axis (arcsec)","Minor Axis Error (arcsec)",'Region Size (pixels)'])

    for file_name in os.listdir(folder_path):
        if file_name.endswith(".txt"):
            file_path = os.path.join(folder_path, file_name)
            process_text_file(file_path, sheet)

    wb.save(output_file)

if __name__ == "__main__":
    input_folder_path = ""
    output_file_path = ""
    process_folder(input_folder_path, output_file_path)
