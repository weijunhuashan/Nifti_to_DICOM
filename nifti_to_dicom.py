#!/usr/bin/env python3
"""Convert cardiac magnitude NIfTI volumes to classic single-frame MR DICOM.
Dependencies: numpy nibabel pydicom>=3.0
Input order for multiple files is EXACT command-line order. No registration/resampling.
Multi-phase export requires timing; --rr-ms assumes uniform phase sampling.
--single-phase [N] exports one phase (1-based, default 1), without requiring timing.
Default output rearranges in-plane pixels so anterior points towards the top.
Viewer display presets may still override the stored orientation.
Not for velocity/phase maps, DWI, segmentation, or Enhanced MR export.
"""
import argparse
import copy
import json
from datetime import datetime
from pathlib import Path
import numpy as np
import nibabel as nib
import pydicom
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import ExplicitVRLittleEndian, MRImageStorage, generate_uid


def dsnum(v):
    return format(float(v), '.10g')


def display_layout(aff, shape, anterior_up=True):
    """Map output (column,row,k) indices to original (i,j,k), without interpolation.

    LPS +Y is posterior. Increasing output row must point posterior so A is up.
    Select the in-plane axis with largest AP projection as the row axis.
    """
    mapping = np.eye(4)
    if anterior_up:
        spacing = np.linalg.norm(aff[:3, :3], axis=0)
        ap = aff[1, :2] / spacing[:2]
        row_axis = int(np.argmax(np.abs(ap)))
        if abs(ap[row_axis]) < 1e-6:
            raise ValueError('Anterior direction is perpendicular to the slice plane; '
                             'A cannot be at the top without reslicing. Use --native-layout '
                             'or resample to an axial/oblique plane first.')
        col_axis = 1 - row_axis
        mapping[:2, :2] = 0
        mapping[col_axis, 0] = 1
        sign = 1 if ap[row_axis] > 0 else -1
        mapping[row_axis, 1] = sign
        if sign < 0:
            mapping[row_axis, 3] = shape[row_axis] - 1
    output_aff = aff @ mapping
    col_axis = int(np.argmax(np.abs(mapping[:2, 0])))
    row_axis = int(np.argmax(np.abs(mapping[:2, 1])))
    output_shape = (shape[col_axis], shape[row_axis], shape[2])
    return mapping, output_aff, output_shape


def arrange_slice(array_ij, mapping):
    col_axis = int(np.argmax(np.abs(mapping[:2, 0])))
    row_axis = int(np.argmax(np.abs(mapping[:2, 1])))
    pixels = np.transpose(array_ij, (row_axis, col_axis))
    if mapping[row_axis, 1] < 0:
        pixels = pixels[::-1, :]
    if mapping[col_axis, 0] < 0:
        pixels = pixels[:, ::-1]
    return pixels


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', nargs='+', required=True, help='One 4D NIfTI OR ordered 3D NIfTI files')
    p.add_argument('--reference', required=True, help='Same-patient reference DICOM')
    p.add_argument('--output', required=True, help='New or empty output directory')
    tm = p.add_mutually_exclusive_group()
    tm.add_argument('--trigger-times-ms', help='Comma-separated delay after R wave, one per phase')
    tm.add_argument('--rr-ms', type=float, help='Nominal RR interval; ASSUMES uniformly spaced phases')
    p.add_argument('--single-phase', nargs='?', const=1, type=int, metavar='N',
                   help='Export one phase only, 1-based (default: first); timing is optional')
    p.add_argument('--start-ms', type=float, default=0, help='First phase delay, used with --rr-ms')
    p.add_argument('--organization', choices=['single-series', 'per-phase'], default='single-series')
    p.add_argument('--description', default='Cardiac NIfTI derived')
    p.add_argument('--series-number', type=int, default=900)
    p.add_argument('--same-frame-of-reference', action='store_true', help='Only if NIfTI retains reference patient coordinates')
    p.add_argument('--assume-mm', action='store_true', help='Explicitly treat unknown NIfTI spatial units as mm')
    p.add_argument('--copy-mr-parameters', action='store_true', help='Copy acquisition parameters ONLY if reference matches source acquisition')
    p.add_argument('--native-layout', action='store_true', help='Keep original in-plane voxel ordering instead of anterior-up')
    a = p.parse_args()
    ref = pydicom.dcmread(a.reference, stop_before_pixels=True)
    out = Path(a.output)
    if out.exists() and any(out.iterdir()):
        p.error('Output directory must be empty; existing output is never overwritten.')
    images = [nib.load(f) for f in a.input]
    if len(images) == 1 and len(images[0].shape) == 4:
        im = images[0]
        frames = [(im, t) for t in range(im.shape[3])]
    elif all(len(im.shape) == 3 for im in images):
        frames = [(im, None) for im in images]
    else:
        p.error('Use one 4D file or a list of 3D files.')
    source_nt = len(frames)
    if source_nt == 0:
        p.error('Input contains no phases.')
    single_phase = a.single_phase is not None
    selected = a.single_phase - 1 if single_phase else None
    if single_phase and not 0 <= selected < source_nt:
        p.error(f'--single-phase must be between 1 and {source_nt}.')
    if not single_phase and source_nt < 2:
        p.error('For a single 3D volume or singleton 4D input, add --single-phase.')
    if not single_phase and a.trigger_times_ms is None and a.rr_ms is None:
        p.error('Multi-phase export requires --rr-ms or --trigger-times-ms.')
    if not np.isfinite(a.start_ms) or a.start_ms < 0:
        p.error('--start-ms must be finite and nonnegative.')
    if a.start_ms != 0 and a.rr_ms is None:
        p.error('--start-ms is only used with --rr-ms.')
    if a.trigger_times_ms is not None:
        try:
            times = np.array([float(x) for x in a.trigger_times_ms.split(',')])
        except ValueError:
            p.error('Trigger times must be comma-separated numbers in milliseconds.')
        # Single export accepts either the selected phase time or the full input list.
        allowed_lengths = {1, source_nt} if single_phase else {source_nt}
        if (len(times) not in allowed_lengths or not np.isfinite(times).all()
                or np.any(times < 0) or np.any(np.diff(times) <= 0)):
            p.error('Trigger times must be finite, nonnegative and increasing; provide '
                    'one per input phase, or one selected-phase time with --single-phase.')
        if single_phase:
            times = times[[selected]] if len(times) == source_nt else times
    elif a.rr_ms is not None:
        if not np.isfinite(a.rr_ms) or a.rr_ms <= 0:
            p.error('--rr-ms must be finite and positive.')
        times = a.start_ms + np.arange(source_nt) * a.rr_ms / source_nt
        if times[-1] >= a.rr_ms:
            p.error('Last phase must be within RR interval; check --start-ms.')
        if single_phase:
            times = times[[selected]]
    else:
        times = None  # Unknown cardiac timing is not invented or copied from reference.
    if single_phase:
        frames = [frames[selected]]
    nt = len(frames)

    # Nibabel affine maps voxel (i,j,k) to RAS; DICOM patient coordinates are LPS.
    def geometry(im):
        if int(im.header['sform_code']) == 0 and int(im.header['qform_code']) == 0:
            raise ValueError('NIfTI has no coded sform/qform; patient geometry is unknown.')
        unit = im.header.get_xyzt_units()[0]
        factors = {'mm': 1., 'meter': 1000., 'micron': .001}
        if unit not in factors and not a.assume_mm:
            raise ValueError('Unknown spatial units: use --assume-mm only after checking.')
        aff = np.diag([-1., -1., 1., 1.]) @ im.affine
        aff[:3, :] *= factors.get(unit, 1.)
        sp = np.linalg.norm(aff[:3, :3], axis=0)
        if not np.isfinite(aff).all() or np.any(sp <= 0):
            raise ValueError('Invalid affine.')
        directions = aff[:3, :3] / sp
        if not np.allclose(directions.T @ directions, np.eye(3), atol=1e-5):
            raise ValueError('Sheared affine: resample to an orthogonal voxel grid first.')
        return aff, sp, directions

    aff, sp, directions = geometry(frames[0][0])
    shape = frames[0][0].shape[:3]
    if max(shape[:2]) > 65535:
        p.error('Rows/Columns exceed classic DICOM limits.')
    mapping, output_aff, output_shape = display_layout(aff, shape, not a.native_layout)
    output_sp = np.linalg.norm(output_aff[:3, :3], axis=0)
    output_dirs = output_aff[:3, :3] / output_sp
    lo, hi = np.inf, -np.inf
    def volume(im, t):
        # ArrayProxy applies NIfTI scl_slope/scl_inter automatically.
        return np.asarray(im.dataobj if t is None else im.dataobj[..., t], dtype=np.float64)
    for im, t in frames:
        af, _, _ = geometry(im)
        if im.shape[:3] != shape or not np.allclose(af, aff, atol=1e-4, rtol=0):
            p.error('All phases must share shape and physical geometry. Register/resample first.')
        v = volume(im, t)
        if not np.isfinite(v).all():
            p.error('NaN/Inf present; clean the data explicitly before conversion.')
        lo, hi = min(lo, float(v.min())), max(hi, float(v.max()))
    # One global linear encoding for ALL phases; reconstructed value = stored*slope+intercept.
    slope = float(dsnum((hi-lo)/65535 if hi > lo else 1.))
    intercept = float(dsnum(lo))
    study = str(getattr(ref, 'StudyInstanceUID', '') or generate_uid())
    frame_uid = generate_uid()
    if a.same_frame_of_reference:
        frame_uid = str(getattr(ref, 'FrameOfReferenceUID', ''))
        if not frame_uid:
            p.error('Reference has no FrameOfReferenceUID.')
    series = [generate_uid() for _ in range(nt if a.organization == 'per-phase' else 1)]
    now = datetime.now()
    copied = ('PatientName PatientID PatientBirthDate PatientSex SpecificCharacterSet '
              'StudyDate StudyTime ReferringPhysicianName StudyID AccessionNumber StudyDescription '
              'PatientPosition').split()
    mr_tags = ('ScanningSequence SequenceVariant ScanOptions MRAcquisitionType RepetitionTime '
               'EchoTime InversionTime EchoTrainLength FlipAngle MagneticFieldStrength '
               'ImagingFrequency ImagedNucleus').split()
    out.mkdir(parents=True, exist_ok=True)
    normal = np.cross(output_dirs[:, 0], output_dirs[:, 1])
    for t, (im, idx) in enumerate(frames):
        v = volume(im, idx)
        folder = out / f'phase_{t+1:03d}'
        folder.mkdir()
        for k in range(shape[2]):
            uid = generate_uid()
            meta = FileMetaDataset()
            meta.MediaStorageSOPClassUID = MRImageStorage
            meta.MediaStorageSOPInstanceUID = uid
            meta.TransferSyntaxUID = ExplicitVRLittleEndian
            filename = folder / f'IM_{k+1:05d}.dcm'
            d = FileDataset(str(filename), {}, file_meta=meta, preamble=b'\0'*128)
            for name in copied + (mr_tags if a.copy_mr_parameters else []):
                if name in ref:
                    d.add(copy.deepcopy(ref.data_element(name)))
            for name in ('PatientName PatientID PatientBirthDate PatientSex StudyDate StudyTime '
                         'ReferringPhysicianName StudyID AccessionNumber ScanningSequence SequenceVariant '
                         'ScanOptions MRAcquisitionType RepetitionTime EchoTime EchoTrainLength').split():
                if name not in d:
                    setattr(d, name, '')
            # Type 1 MR fields: unknown sequence represented by OTHER/NONE.
            if not d.ScanningSequence:
                d.ScanningSequence = 'OTHER'
            if not d.SequenceVariant:
                d.SequenceVariant = 'NONE'
            d.SOPClassUID, d.SOPInstanceUID = MRImageStorage, uid
            d.StudyInstanceUID, d.FrameOfReferenceUID = study, frame_uid
            d.SeriesInstanceUID = series[t if a.organization == 'per-phase' else 0]
            d.Modality = 'MR'
            d.ImageType = ['DERIVED', 'SECONDARY', 'OTHER']
            d.SeriesNumber = a.series_number + (t if a.organization == 'per-phase' else 0)
            d.SeriesDescription = (a.description + (f' phase {t+1:03d}' if a.organization == 'per-phase' else ''))[:64]
            d.InstanceNumber = k+1 if a.organization == 'per-phase' else t*shape[2]+k+1
            d.AcquisitionNumber = t+1
            if nt > 1:
                d.TemporalPositionIdentifier = t+1
                d.NumberOfTemporalPositions = nt
                if np.allclose(np.diff(times), np.diff(times)[0]):
                    d.TemporalResolution = dsnum(times[1]-times[0])
            if times is not None:
                d.TriggerTime = dsnum(times[t])
            if a.rr_ms is not None:
                d.NominalInterval = dsnum(a.rr_ms)
            d.SeriesDate = d.ContentDate = d.InstanceCreationDate = now.strftime('%Y%m%d')
            d.SeriesTime = d.ContentTime = d.InstanceCreationTime = now.strftime('%H%M%S.%f')
            d.DerivationDescription = ('NIfTI conversion; global uint16 linear quantization; '
                                       + ('cardiac timing unknown.' if times is None else
                                          'uniform RR timing estimate.' if a.rr_ms is not None else
                                          'user-supplied trigger times.'))
            d.Manufacturer = 'NIfTI conversion'
            d.PositionReferenceIndicator = ''
            d.Rows, d.Columns = output_shape[1], output_shape[0]
            d.PixelSpacing = [dsnum(output_sp[1]), dsnum(output_sp[0])]
            d.SliceThickness = dsnum(sp[2])  # reconstruction voxel width, not recovered acquisition thickness
            d.SpacingBetweenSlices = dsnum(sp[2])
            pos = (output_aff @ np.array([0., 0., k, 1.]))[:3]
            d.ImagePositionPatient = [dsnum(x) for x in pos]
            d.ImageOrientationPatient = [dsnum(x) for x in np.r_[output_dirs[:, 0], output_dirs[:, 1]]]
            d.SliceLocation = dsnum(pos @ normal)
            d.SamplesPerPixel = 1
            d.PhotometricInterpretation = 'MONOCHROME2'
            d.BitsAllocated = d.BitsStored = 16
            d.HighBit = 15
            d.PixelRepresentation = 0
            d.RescaleSlope, d.RescaleIntercept = dsnum(slope), dsnum(intercept)
            d.RescaleType = 'US'
            d.WindowCenter = dsnum((lo+hi)/2)
            d.WindowWidth = dsnum(max(hi-lo, 1))
            # Rearrange pixels with the EXACT index mapping used for output geometry.
            pixels = np.clip(np.rint((arrange_slice(v[:, :, k], mapping)-intercept)/slope), 0, 65535).astype('<u2')
            d.PixelData = pixels.tobytes()
            d.save_as(filename, enforce_file_format=True)
        timing_text = 'unknown' if times is None else f'{times[t]:g} ms'
        print(f'Phase {t+1}/{nt}: {shape[2]} slices; TriggerTime={timing_text}')
    report = dict(shape=list(shape), phases=nt, files=int(shape[2]*nt),
                  organization=a.organization, trigger_times_ms=None if times is None else times.tolist(),
                  single_phase=single_phase, input_phase_count=source_nt,
                  selected_input_phase=selected+1 if single_phase else None,
                  timing_assumption=('unknown' if times is None else
                                     'uniform RR sampling' if a.rr_ms is not None else 'explicit user times'),
                  affine_lps_mm=aff.tolist(), output_affine_lps_mm=output_aff.tolist(),
                  output_shape=list(output_shape), output_to_input_indices=mapping.tolist(),
                  anterior_up=not a.native_layout, slope=slope, intercept=intercept,
                  max_quantization_error=slope/2, series_uids=series,
                  frame_of_reference_reused=a.same_frame_of_reference)
    (out/'conversion_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('Done. Re-import and verify geometry, intensities, and temporal grouping in target viewer.')


if __name__ == '__main__':
    main()
