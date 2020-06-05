import numpy as np
from scipy.signal.windows import blackman, hann
from scipy.signal import find_peaks
import random
from scipy import interpolate as interp
import matplotlib.pyplot as plt

SILENCE_THRESHOLD = 400.0/32767.0


def silence(x, silence_threshold=SILENCE_THRESHOLD):
    """
    Return boolean indicating if the window is silent (if max of abs(signal) is below the threshold)
    :param x: window signal
    :param silence_threshold: threshold
    :return: boolean, silent (True) or not Silent (False)
    """
    amp_max = np.max(np.abs(x))
    return amp_max < silence_threshold


def build_notes_vector(key, n_oct=4):
    """
    Construct two arrays containing frequency and notes names of notes present in specified key
    :param key: Name of key, has the form 'A' or 'Ab' or 'A#', letters from A to G
    :param n_oct: Number of octave to consider (default 4)
    :return: notes: array containing all frequency of notes in key
            notes_str_ex: array containing all notes names of notes in key
    """
    if key == 'chromatic':
        notes_str = np.asarray(['A', 'A#', 'B', 'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#'], dtype='<U3')
        notes_str_ex = np.tile(notes_str, (1, n_oct + 1))[0]
        octave_index = np.asarray([1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2])
        for idx in range(len(notes_str_ex)):
            incr = idx // 12
            notes_str_ex[idx] = str(notes_str_ex[idx]) + str(octave_index[idx % 12] + incr)
        n_extended = np.arange(len(notes_str_ex))
        notes = np.asarray(55.0 * 2.0 ** (n_extended / 12.0))
        return notes, np.asarray(notes_str_ex)

    else:
        # Name of notes
        if 'b' in key:
            notes_str = ['A', 'Bb', 'B', 'C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab']
            start_idx = notes_str.index(key)
            notes_str = np.asarray(notes_str)
        else:
            notes_str = ['A', 'A#', 'B', 'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#']
            start_idx = notes_str.index(key)
            notes_str = np.asarray(notes_str)

        key_notes = np.sort(np.mod(start_idx + np.array([0, 2, 4, 5, 7, 9, 11]), 12))

    n_extended = np.array([])

    for i in range(n_oct + 1):
        n_extended = np.concatenate((n_extended, 12 * i + key_notes), axis=0)

    all_notes_str_ex = []

    for i in range(n_oct + 1):
        all_notes_str_ex = all_notes_str_ex + [n for n in notes_str]

    for i in range(len(all_notes_str_ex)):
        if i < 3:
            all_notes_str_ex[i] = all_notes_str_ex[i] + str(1)
        else:
            all_notes_str_ex[i] = all_notes_str_ex[i] + str(int((i - 3) / 12 + 2))

    all_notes_str_ex = np.asarray(all_notes_str_ex)
    notes_str_ex = all_notes_str_ex[n_extended.astype(np.int8)]
    notes_str_ex = np.asarray(notes_str_ex)

    # Notes for our table of notes, starting at 55Hz (A1)
    # The factor between each semitone is 2^(1/12)
    notes = np.asarray(55.0 * 2.0 ** (n_extended / 12.0))
    return notes, notes_str_ex


def biased_auto_corr(x, p):
    """
    Compute the biased autocorrelation for computing LPC coefficients
    :param x: frame signal
    :param p: order
    :return: autocorrelation vector
    """
    # compute the biased autocorrelation for x up to lag p
    # Code from LCAV gitbook on LPC
    L = len(x)
    r = np.zeros(p + 1)
    for m in range(0, p + 1):
        for n in range(0, L - m):
            r[m] += x[n] * x[n + m]
        r[m] /= float(L)
    return r


def biased_auto_corr_eff(x, p):
    """
    Compute the biased autocorrelation for computing LPC coefficients efficiently with vectorized operations
    :param x: frame signal
    :param p: order
    :return: autocorrelation vector
    """
    L = len(x)
    r = np.zeros(p + 1)

    for m in range(0, p + 1):
        r[m] = np.sum(x[:L-m] * x[m:])/float(L)


    return r


def levinson_durbin(r, p):
    """
    solve the toeplitz system using the Levinson-Durbin algorithm
    Code from LCAV gitbook on LPC

    :param r: biaised autocorrelation
    :param p: order
    :return: array of coefficient
    """
    if abs(r[0] < 1e-10):
        a = np.zeros(p+1)
        a[0] = 1.0
        return a

    g = r[1] / r[0]
    a = np.array([g])
    v = (1. - g * g) * r[0]

    for i in range(1, p):
        g = (r[i + 1] - np.dot(a, r[1:i + 1])) / v
        a = np.r_[g, a - g * a[i - 1::-1]]
        v *= 1. - g * g
    # return the coefficients of the A(z) filter
    return np.r_[1, -a[::-1]]


def lpc(x, p):
    """
    Compute the LPC coefficients on a given signal frame
    :param x: frame signal (numpy 1D array)
    :param p: order for LPC
    :return: LPC coefficient (numpy 1D array of size p+1)
    """
    # compute p LPC coefficients for a speech segment
    # Code from LCAV gitbook on LPC
    return levinson_durbin(biased_auto_corr_eff(x, p), p)


def yin_algo(x, i, yin_temp, w_len, f_s, f_min, f_max, tol):
    """
    Compute the pitch on the window on x of length yin_len starting at index i using yin algorithm
    :param x: full input signal
    :param i: index
    :param yin_temp: array for algo (to precise)
    :param w_len: length of the window of interest
    :param f_s: sampling frequency
    :param f_min: minimum frequency to detect
    :param f_max: maximum frequency to detect
    :return: pitch in Hz
    """
    tau_max = int(np.round(1 / f_min * f_s))

    # Array to work with (no allocation of mem since no change in size compared to previous iterations)
    x_frame = x[i - tau_max: i + w_len]
    yin_temp = yin_temp * 0

    # Set pitch to 5.0, will be modified if a pitch is detected
    pitch = 5.0

    for tau in range(tau_max):
        yin_temp[tau] = np.sum(((x_frame[0:w_len] - x_frame[tau: w_len + tau]) ** 2))

    tmp = 0
    yin_temp[0] = 1

    for tau in range(1, tau_max):
        tmp = tmp + yin_temp[tau]
        if tmp == 0:
            return pitch
        yin_temp[tau] = yin_temp[tau] * tau / tmp

    tau = int(f_s / f_max)

    while tau < tau_max:
        if yin_temp[tau] < tol:
            while yin_temp[tau + 1] < yin_temp[tau]:
                tau = tau + 1
                if tau + 1 >= tau_max:
                    break

            pitch = f_s / tau

            break
        else:
            tau = tau + 1

    return pitch


def pitch_marks(x_frame, pitch, prev_marks, prev_pitch, prev_voiced_pitch, w_len, hop, f_s, delta, valley=True):
    """
    Computes the pitch marks position given the pitch and the previous pitch marks positions
    :param x_frame: Current window on which to find analysis pitch marks
    :param pitch: current pitch (detected on the current window)
    :param prev_marks: previous pitch marks array
    :param prev_pitch: previous pitch (detected on previous window)
    :param prev_voiced_pitch: last voiced pitch detected (last value of pitch detected (not necessarily on prev window))
    :param w_len: length of window
    :param hop: hop size
    :param f_s: sampling frequency
    :param delta: tolerance number to search for extremum
    :param valley: (bool) if True, put pitch marks on minima, otherwise on maxima
    :return: array with indices of pitch marks
    """

    def arg_ext(frame, valley):
        if valley:
            return np.argmin(frame)
        else:
            return np.argmax(frame)

    # Pitch period (in number of samples)
    T_prev = int(f_s / prev_voiced_pitch)

    # Set prev_marks to curr referential of idx
    prev_marks = prev_marks - hop
    n_marks_ov = np.sum(prev_marks >= 0)

    search_left = False

    # If window is voiced (if a pitch is detected)
    if pitch > 10:
        T = int(f_s / pitch)

        # Close search limit
        sw_c = int(np.floor(delta * T))
        # Far search limit
        sw_f = int(np.ceil((2.0 - delta) * T))

        # Placement of the first mark

        # If last window was voiced
        if prev_pitch > 10:
            # No need to search left after since we are sure it is the fist mark of the window

            # If no previous marks in overlapping parts
            if n_marks_ov == 0:
                # Assure to stay in bounds
                last_mark = prev_marks[-1]
                l_lim = np.max((last_mark + int(np.min((sw_c, np.floor(delta * T_prev)))), 0))
                r_lim = np.min((last_mark + int(np.max((sw_f, np.ceil((2.0 - delta) * T_prev)))), w_len))

                t = arg_ext(x_frame[l_lim: r_lim], valley)
            else:
                # Take first of pitch mark in overlap as first t
                t = prev_marks[-n_marks_ov]
        else:
            # If no info before (previous window was unvoiced): discard any info about previous pitch mark
            # Find arg_ext on all window, and then search right and left
            search_left = True
            t = arg_ext(x_frame, valley)

        marks = np.array([t])

        # Search to the right
        while marks[-1] + sw_c < w_len:
            # Check if all search area is contained in frame
            if marks[-1] + sw_f < w_len:
                extr_idx = arg_ext(x_frame[marks[-1] + sw_c: marks[-1] + sw_f], valley)
                marks = np.append(marks, marks[-1] + sw_c + extr_idx)
            else:
                if marks[-1] + T < w_len:
                    extr_idx = arg_ext(x_frame[marks[-1] + sw_c:], valley)
                    marks = np.append(marks, marks[-1] + sw_c + extr_idx)
                    break
                else:
                    break

        # Search to the left
        if search_left:
            while marks[0] - sw_c > 0:
                # Check if all search area is contained in frame
                if marks[0] - sw_f > 0:
                    extr_idx = arg_ext(x_frame[marks[0] - sw_f: marks[0] - sw_c], valley)
                    marks = np.insert(marks, 0, marks[0] - sw_f + extr_idx)
                else:
                    if marks[0] - T >= 0:
                        extr_idx = arg_ext(x_frame[0: marks[0] - sw_c], valley)
                        marks = np.insert(marks, 0, extr_idx)
                        break
                    else:
                        break

    # if window is not voiced (no pitch detected)
    else:
        # If prev_marks empty and window not voiced
        if prev_marks.size == 0:
            return np.array([])
        # Put pitch mark at a constant rate w.r.t to prev_voiced_pitch
        # Place first mark(s)
        if n_marks_ov > 0:
            marks = prev_marks[prev_marks >= 0]
        else:
            marks = np.array([prev_marks[-1] + T_prev])

        while marks[-1] + T_prev < w_len:
            marks = np.append(marks, marks[-1] + T_prev)

    return marks


def synthesis_pitch_marks(pitch, prev_pitch, prev_voiced_pitch, an_marks, prev_st_marks, beta,
                     w_len, hop, f_s):
    if pitch > 10:
        pitch_new = beta * pitch
        T = int(f_s / pitch)

    else:
        pitch_new = beta * prev_voiced_pitch
        T = int(f_s/prev_voiced_pitch)

    T_new = int(f_s / pitch_new)

    # Set idx to same referential
    prev_st_marks = prev_st_marks - hop
    n_marks_ov = np.sum(prev_st_marks >= 0)

    N = w_len // T_new

    # if current window is voiced
    if pitch > 10:
        # if previous is also voiced
        if prev_pitch > 10:
            # Maintain continuity of synthesis marks if curr is voiced and prev is voiced
            if n_marks_ov > 0:
                st_marks = prev_st_marks[-n_marks_ov] + np.arange(0, N + 2) * T_new
            else:
                st_marks = prev_st_marks[-1] + np.arange(1, N + 2) * T_new

        # if previous unvoiced
        else:
            st_marks = an_marks[0] + np.arange(-3, N + 2) * T_new
    else:

        if prev_st_marks.size == 0:
            return np.array([])

        # if previous unvoiced also
        if prev_pitch < 10:
            # Maintain continuity of synthesis marks if curr is unvoiced and prev is unvoiced
            if n_marks_ov > 0:
                st_marks = prev_st_marks[-n_marks_ov] + np.arange(0, N + 2) * T_new
            else:
                st_marks = prev_st_marks[-1] + np.arange(1, N + 2) * T_new

        # if previous is voiced
        else:
            # Maintain continuity of synthesis marks if curr is unvoiced and prev is voiced
            if n_marks_ov > 0:
                st_marks = prev_st_marks[-n_marks_ov] + np.arange(0, N + 2) * T_new
            else:
                st_marks = prev_st_marks[-1] + np.arange(1, N + 2) * T_new

    # Remove marks outside of frame
    st_marks = st_marks[np.logical_and(st_marks >= 0, st_marks < w_len)]

    return st_marks


def pitch_shift(e, out_window, pitch, prev_voiced_pitch, an_marks, st_marks, beta,
                     w_len, f_s, tau_max):
    """
    Place synthesis marks and apply psola on residual signal
    :param e: residual frame
    :param out_window: residual window
    :param pitch: current pitch
    :param prev_voiced_pitch: previous voiced pitch
    :param an_marks: array of analysis pitch marks
    :param st_marks: synthesis pitch marks
    :param beta: shift factor
    :param w_len: frame length
    :param f_s: sampling frequency
    :param tau_max: max period fs/fmin
    :return: out_window and synthesis marks
    """

    if pitch > 10:
        pitch_new = beta * pitch
        T = int(f_s / pitch)

    else:
        pitch_new = beta * prev_voiced_pitch
        T = int(f_s/prev_voiced_pitch)

    T_new = int(f_s / pitch_new)

    window = hann(2*T+1)

    if st_marks.size == 0:
        return e[tau_max:]

    # for every synthesis mark
    for idx, mark in enumerate(st_marks):
        # Closest analysis marks
        cl_mark = an_marks[np.argmin(np.abs(mark - an_marks))]

        # PSOLA with all the different possible cases

        # If only one pitch mark in frame
        if idx == 0 and st_marks.shape[0] == 1:
            print("Only one st_marks")
            # No window applied
            x_p = e[tau_max + cl_mark - T: np.min((tau_max + cl_mark + T + 1, e.shape[0]))]
            t = np.arange(-T, - T + x_p.shape[0])

            int_func = interp.interp1d(mark + t / beta, x_p, kind='linear', bounds_error=False,
                                       fill_value=0.0, assume_sorted=True)

            out_window[np.max((0, mark - T_new)): np.min((w_len, mark + T_new + 1))] += int_func(np.arange(
                np.max((0, mark - T_new)),
                np.min((w_len, mark + T_new + 1))
            ))

        elif idx == 0:
            if cl_mark + T + 1 < w_len:
                l_window = np.concatenate((np.ones(T), window[T:]))
            else:
                l_window = np.concatenate((np.ones(T), window[T:]))

            x_p = e[tau_max + cl_mark - T: tau_max + cl_mark + T + 1] * l_window

            t = np.arange(-T, T + 1)
            int_func = interp.interp1d(mark + t / beta, x_p, kind='linear', bounds_error=False,
                                       fill_value=0.0, assume_sorted=True)

            out_window[np.max((0, mark - T_new)): mark + T_new + 1] += int_func(
                np.arange(np.max((0, mark - T_new)), mark + T_new + 1))

        elif idx == st_marks.shape[0] - 1:
            x_p = e[tau_max + cl_mark - T: np.min((tau_max + cl_mark + T + 1, e.shape[0]))]

            l = x_p.shape[0]
            lr = l - (T + 1)

            r_window = np.concatenate((window[:T + 1], np.ones(lr)))
            x_p *= r_window

            t = np.arange(-T, lr + 1)
            int_func = interp.interp1d(mark + t / beta, x_p, kind='linear', bounds_error=False,
                                       fill_value=0.0, assume_sorted=True)

            out_window[mark - T_new: np.min((w_len, mark + T_new + 1))] += \
                int_func(np.arange(mark - T_new, np.min((mark + T_new + 1, w_len))))
        else:
            # Extract 2 period centered at cl mark
            e_2 = e[tau_max + cl_mark - T: np.min((tau_max + cl_mark + T + 1, e.shape[0]))]

            if e_2.shape[0] != 2 * T + 1 and e_2.shape[0] != 0:
                window_2 = window[:e_2.shape[0]]
                x_p = e_2 * window_2
                t = np.arange(-T, e_2.shape[0] - T)

            else:
                x_p = e_2 * window
                t = np.arange(-T, T + 1)

            int_func = interp.interp1d(mark + t / beta, x_p, kind='linear', bounds_error=False,
                                       fill_value=0.0, assume_sorted=True)

            out_window[mark - T_new: mark + T_new + 1] += int_func(np.arange(mark - T_new, mark + T_new + 1))

    return out_window


def create_window(w_size, overlap=0.5, type='sine'):
    """
    Return array of size w_size with window function evaluated at each index
    :param w_size: size of window
    :param overlap: 0, 0.5 or 0.75 of overlap
    :param type: type of the windows (sine, hann, hamming, or rect)
    :return: array containing window function evaluated between 0 and w_size-1
    """
    if overlap==0.75:
        overlap_factor = 1.0/np.sqrt(2)
    elif overlap==0.5:
        overlap_factor = 1.0
    elif overlap==0.0:
        w = np.ones(w_size)
        return w
    else:
        raise ValueError('Not valid overlap, should be 0, 0.5 or 0.75')

    n = np.arange(w_size)

    if type == 'sine':
        w = overlap_factor * np.sin((n+0.5)*np.pi/w_size)

    elif type == 'hann':
        w = overlap_factor * np.sin((n+0.5)*np.pi/w_size)**2

    elif type == 'rect':
        w = overlap_factor * np.ones(w_size)

    elif type == 'hamm':
        w = overlap_factor * np.hamming(w_size)

    else:
        raise ValueError('Not valid window type')

    return w