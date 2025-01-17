file_name = "/Users/shivaniramesh/Desktop/17mix_test2.mzxml.gz"  

import sys
from base64 import b64decode
from array import array
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import gzip


def get_spectrum_from_scan(file_name, scan_number):

    xmldoc = gzip.open(file_name)

    mzs, ints = None, None
    scan_found = False
    ns = '{http://sashimi.sourceforge.net/schema/}'
    # Iterate through the file using iterparse
    for event,ele in ET.iterparse(xmldoc):
        if not ns:
            p = ele.tag.find('}')
            if p >= 0:
                ns = ele.tag[:(p+1)]

        # Check for the scan element
        if event == 'end' and ele.tag == ns + 'scan' and int(ele.attrib.get('num', -1)) == scan_number:

            scan_found = True  # Mark the desired scan as found

        # Check for the peaks element within the matched scan
        if event == 'end' and ele.tag == ns + 'peaks' and scan_found:

            if ele.text is None or not ele.text.strip():
                print(f"No peaks found for scan {scan_number}.")
                return None, None
            
            # Decode peak data
            peaks = array("f", b64decode(ele.text.strip()))
            if sys.byteorder != "big":
                peaks.byteswap()
            mzs = peaks[::2]  # m/z values (even indices)

            ints = peaks[1::2]  # Intensity values (odd indices)

            scan_found = False  # Reset for next scans if any
            ele.clear()  # Free memory for parsed element
            break  # Exit once we find the desired scan
        
        ele.clear()  # Free memory for other elements

    if mzs is None:
        print(f"Scan {scan_number} not found in the mzXML file.")
    return mzs, ints


def compute_ion_masses(peptide):

    aa_mass = {
        'A': 71.04, 'R': 156.10, 'N': 114.04, 'D': 115.03,
        'C': 103.01, 'E': 129.04, 'Q': 128.06, 'G': 57.02,
        'H': 137.06, 'I': 113.08, 'L': 113.08, 'K': 128.10,
        'M': 131.04, 'F': 147.07, 'P': 97.05, 'S': 87.03,
        'T': 101.05, 'W': 186.08, 'Y': 163.06, 'V': 99.07
    }


    #initialize aa ion mass list
    b_ions = []
    y_ions = []

    b_mass = 1 #proton gained at cleavage of N-term
    y_mass = 1 + 18  #Water molecule for y-ions at C-term (COO- +H + H2O )
    
    for i, aa in enumerate(peptide):
        
        b_mass += aa_mass[aa]
        b_ions.append(b_mass)  
        
        y_mass += aa_mass[peptide[-(i + 1)]] #backwards
        y_ions.append(y_mass)

        if not b_ions or not y_ions:
            print("No valid amino acids found in peptide sequence.")

    return b_ions, y_ions


def plot_peptide_match(mzs, ints, b_ions, y_ions, peptide, scan_number, tolerance= 0.5):
    #normalize to % of max ints
    max_intensity = max(ints)
    #stem plot to show experimental spectrum
    plt.figure(figsize=(14, 8))
    plt.stem(mzs, ints, linefmt="green", markerfmt=" ", basefmt=" ", label="Experimental Spectrum")
   


    def find_closest_peak(theoretical_mz, mzs, ints, tolerance):
        #filter for absolute difference of m/z <=  tolerance AND > thershold : save as (mz, int) coordinate
        filtered_peaks = [
            (mz, intensity) for mz, intensity in zip(mzs, ints)
            if abs(mz - theoretical_mz) <= tolerance and intensity > 5  # Threshold: 5% of max intensity
        ]
        return min(filtered_peaks, key=lambda x: abs(x[0] - theoretical_mz), default=None) #return closest peak, or None
    
    # Annotate b-ion peaks
    for i, b_mz in enumerate(b_ions):
        closest_peak = find_closest_peak(b_mz, mzs, ints, tolerance)  # Find the closest peak
        if closest_peak:
            closest_mz, closest_intensity = closest_peak 

            #annotate ion peak and add text , print matches in console
            plt.text(closest_mz, closest_intensity + 2, f"b{i+1}", color="blue", fontsize=10)
            plt.stem([closest_mz], [closest_intensity], linefmt="blue", markerfmt=" ", basefmt=" ")
            print("MATCH -- [b", i+1, "] : b_mz", b_mz , " = ", closest_mz)
    
    # Annotate y-ion peaks
    for i, y_mz in enumerate(y_ions):
        closest_peak = find_closest_peak(y_mz, mzs, ints, tolerance)
        if closest_peak:
            closest_mz, closest_intensity = closest_peak

            plt.text(closest_mz, closest_intensity + 4, f"y{i+1}", color="red", fontsize=10)
            plt.stem([closest_mz], [closest_intensity], linefmt="red", markerfmt=" ", basefmt=" ")

            print("MATCH -- [y", i+1, "] : y_mz", y_mz , " = ", closest_mz)


    plt.title(f"Annotated Peptide Fragmentation Spectrum (Scan {scan_number})\nPeptide: {peptide}", fontsize=14)
    plt.xlabel("m/z (Mass-to-Charge Ratio)", fontsize=12)
    plt.ylabel("Relative Abundance (%)", fontsize=12)
    plt.xlim(min(mzs) - 50, max(mzs) + 50)
    plt.ylim(0, max_intensity + 50)  
    plt.legend(fontsize=10, loc="upper right")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


def main():

    file_name = sys.argv[1]

    scan_number = int(sys.argv[2])
    peptide = sys.argv[3]
    
    # Compute b-ion and y-ion m/z values
    try:
        b_ions, y_ions = compute_ion_masses(peptide)
    except KeyError:
        print("Error: Invalid amino acid in the peptide sequence.")
        return
    
    # Extract spectrum 
    mzs, ints = get_spectrum_from_scan(file_name, scan_number)
    if mzs is None:
        print("Scan ", scan_number, "not found in ",file_name, ".")
        return
    
    # Plot
    plot_peptide_match(mzs, ints, b_ions, y_ions, peptide, scan_number)

if __name__ == "__main__":
    main()
