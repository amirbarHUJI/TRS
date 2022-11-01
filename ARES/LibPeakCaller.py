import logging
import numpy as np

from scipy.signal import find_peaks


from ARES.Globals import Strand, Peak


class LibPeakCaller(object):

    def get_signal_dct(self, signal_counts, normalize_func=None):
        signal_dct = {}

        for chr_name, chr_len in signal_counts.Chromosomes:

            signal_dct[chr_name] = {
                Strand.FWD: [],
                Strand.REV: []
            }

            for strand in (Strand.FWD, Strand.REV):
                signal_dct[chr_name][strand] = \
                    signal_counts[chr_name, strand, 1, chr_len]

                if normalize_func:
                    signal_dct[chr_name][strand] = \
                        normalize_func(signal_dct[chr_name][strand])

        return signal_dct

    def get_peak_dct(
        self,
        signal_dct,
        prominence,
        width,
        rel_height,
        distance,
        min_height,
    ):

        peak_dct = {}

        for chr_name in signal_dct:
            for strand, signal_array in signal_dct[chr_name].items():
                peak_dct[(chr_name, strand.value)] = self.call_strand_peaks(
                    signal_array=signal_array,
                    chr_name=chr_name,
                    strand=strand,
                    prominence=prominence,
                    width=width,
                    rel_height=rel_height,
                    distance=distance,
                    min_height=min_height)

        return peak_dct

    def call_strand_peaks(
        self,
        signal_array,
        chr_name,
        strand,
        prominence,
        width,
        rel_height,
        distance,
        min_height,
    ):
        peaks, properties = find_peaks(
            height=min_height,
            x=signal_array,
            prominence=prominence,
            width=width,
            rel_height=rel_height,
            distance=distance
        )

        starts = np.arange(0, len(signal_array))
        left_ips = properties["left_ips"]
        right_ips = properties["right_ips"]
        prominences = properties["prominences"]

        xp = np.arange(len(starts))
        ileft_ips = np.interp(left_ips, xp, starts).round().astype(int)
        iright_ips = np.interp(right_ips, xp, starts).round().astype(int)
        ipeaks = np.interp(peaks, xp, starts).round().astype(int)

        idx = ileft_ips <= iright_ips
        ileft_ips = ileft_ips[idx]
        iright_ips = iright_ips[idx]
        ipeaks = ipeaks[idx]
        signals = signal_array[peaks[idx]]

        removed_peaks = sum(~idx)

        if removed_peaks:
            logging.info("Removed %d invalid peaks. " % removed_peaks)

        peak_list = [
            Peak(
                chr_name,
                strand.value,
                start,
                end,
                signal,
                peak
            )
            for start, end, peak, signal in zip(
                ileft_ips, iright_ips, ipeaks, signals)
        ]

        return peak_list
