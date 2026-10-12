import numpy as np

import logging

from astropy.nddata import CCDData, VarianceUncertainty
from astropy import wcs
import astropy.units as u
from astropy.table import Table

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.geom as geom

from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
import lsst.pipe.base.connectionTypes as connTypes

from lsst.pex.exceptions import InvalidParameterError
from lsst.afw.detection import InvalidPsfError



class SMPTaskConnections(PipelineTaskConnections,
                        dimensions=("instrument", "visit", "detector")):
    
    visitImage = connTypes.Input(
        doc="Visit image (calibrated exposure)",
        name="visit_image",
        dimensions=("instrument","visit","detector"),
        storageClass="ExposureF",
        multiple=True,
    )
    
    visitSummary = connTypes.Input(
        doc="Summary of the visits.",
        name="visit_summary",
        dimensions=("instrument", "visit"),
        storageClass="ExposureCatalog",
    )

    scene_modeling_inputs = connTypes.Output(
        doc="Inputs prepared for scene modeling photometry.",
        name="sceneModelingInputs",
        dimensions=("instrument","visit","detector"),
        storageClass="ArrowAstropy",
    )



def _validate_dec(item):
    "Ensure a valid declination value."
    return (item > -90) and (item < 90)

def _validate_ra(item):
    "Ensure a valid right ascention value."
    return (item > 0) and (item < 360)



class SMPTaskConfig(
        PipelineTaskConfig,
        pipelineConnections=SMPTaskConnections):
    """Configuration for CreateTask, prepares scene modeling inputs.
    """

    ra = pexConfig.ListField(
        doc="RA of the SN candidates (degrees).",
        dtype=float,
        maxLength=50,
        itemCheck=_validate_ra,
    )

    dec = pexConfig.ListField(
        doc="Dec of the SN candidates (degrees).",
        dtype=float,
        maxLength=50,
        itemCheck=_validate_dec,
    )
    
    zero_flux_min_mjd = pexConfig.ListField(
        doc="MJD before which the flux at the SN position is zero (OFF visits).",
        dtype=float,
        default=(60310,),
    )

    zero_flux_max_mjd = pexConfig.ListField(
        doc="MJD after which the flux at the SN position is zero (OFF visits).",
        dtype=float,
        default=(65000,),
    )

    cutout_half_size_pixels = pexConfig.Field(
        doc="Half-size of the cutout in pixels (full cutout = 2*N+1 x 2*N+1).",
        dtype=int,
        default=25,
    )


class SMPTask(PipelineTask):
    """Prepare inputs for scene modeling photometry of SNIa candidates.
    """

    _DefaultName = "createSceneModelingInputs"
    ConfigClass = SMPTaskConfig

    def __init__(self, **kwargs):
        super().__init__(**kwargs)   # initialise the mother class PipelineTask
        self.log = logging.getLogger(__name__)   # initialise a Python logger for emitting follow-up messages (info, warning, debug)
        self._scale = 1.0

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, visitImage_refs, visitSummary=None):
        """Execute the pipeline task to prepare inputs for scene modeling photometry.

        Iterates over input visit exposure references, checks sky coverage for each SN 
        candidate, extracts image cutouts centered on the source, computes the 
        local PSF, retrieves photometric calibration ratios, and evaluates the 
        ON/OFF state based on observation epochs (MJD).

        Parameters
        ----------
        visitImage_refs : 'list' of 'lsst.daf.butler._deferredDatasetHandle.DeferredDatasetHandle'
            References to the calibrated visit images ('visit_image') to process.
        visitSummary : 'lsst.afw.table.ExposureCatalog', optional
            Summary catalog containing information about the visits.

        Returns
        -------
        'lsst.pipe.base.Struct'
            A structure containing a dictionary under the key 'sceneModelingInputs' with:
            - 'cutouts': 'list' of 'astropy.nddata.CCDData' objects containing image cutouts, WCS, variances and masks.
            - 'positions': 'list' of 'tuple' (ra, dec) coordinates for each processed target.
            - 'on_off': 'list' of 'bool' flags indicating if the SN is expected to be present (ON) or absent (OFF).
            - 'photometric_ratios': 'list' of 'float' local calibration factors.
        """

        cutouts = []
        positions = []
        on_off_flags = []
        phot_ratios = []
        
        for visit_ref in visitImage_refs:
            visitImage = visit_ref.get()
            wcs = visitImage.getWcs()
            #psf = visitImage.getPsf()
            photoCalib = visitImage.getPhotoCalib()
            bbox = visitImage.getBBox()
            visitInfo = visitImage.getInfo().getVisitInfo()
            mjd = visitInfo.getDate().get()  # MJD as float

            #print(photoCalib)
            #print(f"calibrationMean = {photoCalib.getCalibrationMean()}")
            #print(f"calibrationErr  = {photoCalib.getCalibrationErr()}")
    
            # Determine ON/OFF for this visit
            on_flag = self._compute_on_off(mjd)
    
            for ra, dec in zip(self.config.ra, self.config.dec):
                sky_coord = geom.SpherePoint(ra, dec, geom.degrees)  # convert SN candidate coords in lsst.geom object
                pixel_coord = wcs.skyToPixel(sky_coord)
                pixel_point = geom.Point2I(pixel_coord)
                
                if visitImage.containsSkyCoords(ra*u.degree, dec*u.degree) == False :
                    continue
    
                # 1- Data cutout
                cutout = self._extract_cutout(visitImage, sky_coord)
                if cutout is None:
                    self.log.warning(
                        f"Could not extract cutout for SN at (ra={ra}, dec={dec}).")
                    continue
                ccdData, raw_array = cutout
    
                # 2- PSF at SN position
                cutout_psf = self._compute_psf(visitImage, pixel_coord)  # now done in _extract_cutout
                ccdData.psf = cutout_psf
    
                # 3- PSF derivative
    
                # 4- Photometric ratio (Ri)
                phot_ratio = self._compute_phot_ratio(photoCalib, pixel_coord)
                
                # Append results
                cutouts.append(ccdData)
                positions.append((ra, dec))
                on_off_flags.append(on_flag)
                #psfs.append(cutout_psf)
                phot_ratios.append(phot_ratio)
    
        result_dict = {
            "cutouts": cutouts,
            "positions": positions,
            "on_off": on_off_flags,
            "photometric_ratios": phot_ratios
        }

        return Struct(sceneModelingInputs=result_dict)

    

    ########### Private methods ####################################################################

    def _compute_on_off(self, mjd):
        """ Return True if the SN is expected to be present at this MJD.
        A visit is ON if mjd is between the (min, max) pair in the config.
        It is OFF otherwise (before the explosion or after the SN has faded).
        """
        for mjd_min, mjd_max in zip(self.config.zero_flux_min_mjd, self.config.zero_flux_max_mjd):
            if (mjd_min <= mjd) and (mjd <= mjd_max):
                return True
        return False

    
    def _extract_cutout(self, visitImage, sky_coord):
        """Extract a square cutout array centered on the given sky_coord (SpherePoint object).

        Projects the target sky coordinates to pixel coordinates on the provided 
        visit image, defines a square bounding box using the configured half-size, 
        and extracts the corresponding cutout exposure. It applies photometric 
        calibration, builds a local Astropy WCS approximation, and wraps the 
        resulting data into an Astropy CCDData object.

        Parameters
        ----------
        visitImage : 'lsst.afw.image.Exposure'
            The calibrated visit image exposure from which to extract the cutout.
        sky_coord : 'lsst.geom.SpherePoint'
            The (RA, Dec) coordinates representing the center of the cutout.

        Returns
        -------
        'tuple' of ('astropy.nddata.CCDData', 'numpy.ndarray') or 'None'
            A tuple containing:
            - `ccdData`: The calibrated cutout wrapped as an Astropy CCDData object 
              (including data array, variance uncertainty, mask, local WCS, and metadata).
            - `raw_array`: A copy of the uncalibrated image data array.
            Returns `None` if the calculated cutout box extends outside the valid 
            bounding box of the visit image.
        """
        visit_wcs = visitImage.getWcs()
        pixel_point = visit_wcs.skyToPixel(sky_coord)
        
        x_center = int(round(pixel_point.x))
        y_center = int(round(pixel_point.y))
        
        half_size = self.config.cutout_half_size_pixels
        extent = geom.Extent2I(2 * half_size + 1, 2 * half_size + 1)
        
        min_point = geom.Point2I(x_center - half_size, y_center - half_size)
        cutout_bbox = geom.Box2I(min_point, extent)

        if not visitImage.getBBox().contains(cutout_bbox):
            self.log.warning(
                f"Cutout {cutout_bbox} extends outside the image bounding box {visitImage.getBBox()}."
            )
            return None

        # Extract cutout
        cutout_exposure = visitImage[cutout_bbox]

        # Get photometric calibration info
        photoCalib = visitImage.getPhotoCalib()
        calibCutout = photoCalib.calibrateImage(cutout_exposure.getMaskedImage())

        cutOutMinX = cutout_exposure.getBBox().minX - 1
        cutOutMinY = cutout_exposure.getBBox().minY - 1

        cutoutWcs = wcs.WCS(naxis=2)
        cutoutWcs.array_shape = (cutout_exposure.getBBox().getHeight(),
                                 cutout_exposure.getBBox().getWidth())
        cutoutWcs.wcs.crpix = [pixel_point.x - cutOutMinX, pixel_point.y - cutOutMinY]
        cutoutWcs.wcs.crval = [sky_coord.getRa().asDegrees(),
                               sky_coord.getDec().asDegrees()]
        cutoutWcs.wcs.cd = self._make_local_transform_matrix(visit_wcs, pixel_point, sky_coord)

        #cutout2D = np.array([calibCutout.getImage().array, cutout_exposure.getImage().array])

        # Build a CCDData object with everything
        ccdData = CCDData(
            data=calibCutout.getImage().array,
            uncertainty=VarianceUncertainty(calibCutout.getVariance().array),
            flags=calibCutout.getMask().array,
            wcs=cutoutWcs,
            meta={"cutMinX": cutOutMinX,
                  "cutMinY": cutOutMinY},
            unit=u.nJy)

        return ccdData, cutout_exposure.getImage().array.copy()

        
    def _make_local_transform_matrix(self, wcs, center, skyCenter):
        """Create a local, linear approximation of the wcs transformation matrix.
        (from https://github.com/lsst/ap_association/blob/tickets/DM-52481/python/lsst/ap/association/packageAlerts.py)

        The approximation is created as if the center is at RA=0, DEC=0. All
        comparing x,y coordinate are relative to the position of center. Matrix
        is initially calculated with units arcseconds and then converted to
        degrees. This yields higher precision results due to quirks in AST.

        Parameters
        ----------
        wcs : `lsst.afw.geom.SkyWcs`
            Wcs to approximate
        center : `lsst.geom.Point2D`
            Point at which to evaluate the LocalWcs.
        skyCenter : `lsst.geom.SpherePoint`
            Point on sky to approximate the Wcs.

        Returns
        -------
        localMatrix : `numpy.ndarray`
            Matrix representation the local wcs approximation with units
            degrees.
        """
        blankCDMatrix = [[self._scale, 0], [0, self._scale]]
        localGnomonicWcs = afwGeom.makeSkyWcs(
            center, skyCenter, blankCDMatrix)
        measurementToLocalGnomonic = wcs.getTransform().then(
            localGnomonicWcs.getTransform().inverted()
        )
        localMatrix = measurementToLocalGnomonic.getJacobian(center)
        return localMatrix / 3600
        

    def _compute_psf(self, exposure, pixelCenter, srcId=None):
        """Compute the PSF at a location and catch errors.
        (from https://github.com/lsst/ap_association/blob/tickets/DM-52481/python/lsst/ap/association/packageAlerts.py)

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The image to compute the PSF for.
        pixelCenter : `lsst.geom.Point2D`
            The location on the image to compute the PSF.
        srcId : `int`, optional
            Unique id of DiaSource. Used for when an error occurs extracting
            a cutout.

        Returns
        -------
        cutoutPsf : `numpy.array`
            Array of the PSF values.
        """
        try:
            # use exposure.psf.computeKernelImage to provide PSF centered in the array
            cutoutPsf = exposure.psf.computeKernelImage(pixelCenter).array
        except InvalidParameterError:
            if srcId is not None:
                msg = "Could not calculate PSF for DiaSource with "\
                      "id=%i. InvalidParameterError encountered. Exiting."\
                      % srcId
            else:
                msg = "Could not calculate average PSF for the image"
            self.log.warning(msg)
            cutoutPsf = None
        except InvalidPsfError:
            if srcId is not None:
                msg = "Could not calculate PSF for DiaSource with "\
                      "id=%i. InvalidPsfError encountered. Exiting."\
                      % srcId
            else:
                msg = "Could not calculate average PSF for the image"
            self.log.warning(msg)
            cutoutPsf = None
        # Cast the PSF to float32 to reduce the size of the alert.
        if cutoutPsf is not None:
            cutoutPsf = cutoutPsf.astype(np.float32)

        return cutoutPsf

        
    #def _compute_psf_derivative(self):
    #    return

    
    def _compute_phot_ratio(self, photoCalib, pixel_coord):
        """Return the local photometric calibration factor R_i [nJy / electron] defined as :
            F_physical = R_i * f_instrumental
        
        Parameters
        ----------
        photoCalib : `lsst.afw.image.PhotoCalib`
            Photometric calibration object attached to the visit image.
        pixel_coord : `lsst.geom.Point2D`
            Pixel coordinates of the SN on this detector.

        Returns
        -------
        phot_ratio : `float`
            Local calibration factor in nJy/electron
            (Note: if the mean equals 1.0 with an error of about 5e-4, it means the visitImage
            has already been pre-normalized during dataset preparation, meaning that upstream
            processing has already converted raw instrumental counts into physical units (nJy)
            directly within the pixels, so the calibration factor acts as a neutral identity).
        """
        return float(photoCalib.getLocalCalibration(pixel_coord))