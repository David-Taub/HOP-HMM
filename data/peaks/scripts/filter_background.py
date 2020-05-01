import sys
bed_filepath = sys.argv[1]
peak_length = int(sys.argv[2])
print(peak_length)
with open(bed_filepath, 'r') as f:
    peaks = f.readlines()
peaks = [peak[:-1].split('\t') for peak in peaks]
print(len(peaks))
peaks = [peak for peak in peaks if (int(peak[2]) - int(peak[1])) > peak_length]
print(len(peaks))
peaks = [[peak[0], str(int((int(peak[2]) + int(peak[1]) - peak_length) / 2)),
          str(int((int(peak[2]) + int(peak[1]) + peak_length) / 2))] for peak in peaks]
peaks = ['\t'.join(peak) for peak in peaks]
with open(bed_filepath, 'w') as f:
    f.write('\n'.join(peaks))
