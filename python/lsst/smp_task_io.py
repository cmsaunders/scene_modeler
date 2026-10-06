import numpy as np

import lsst.pex.config as pexConfig
import lsst.afw.image as afwImage
import lsst.geom as geom

from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
import lsst.pipe.base.connectionTypes as connTypes



class CreateConnections(PipelineTaskConnections,
                        dimensions=("instrument", "visit", "detector")):
    
    visitImage = connTypes.Input(
        doc="Visit image (calibrated exposure)",
        name="visit_image",
        dimensions=("instrument","visit","detector"),
        storageClass="ExposureF",
        multiple_inputs=True,
    )
    
    visitSummary = connTypes.Input(
        doc="Summary of the visits.",
        name="visit_summary",
        dimensions=("instrument", "visit"),
        storageClass="ExposureCatalog",
    )

    ########
    scene_modeling_inputs = connTypes.Output(
        doc="Inputs prepared for scene modeling photometry.",
        name="",
        dimensions="",
        storageClass="ArrowAstropy",
    )
    ########


def _validate_dec(item):
    "Ensure a valid declination value."
    return (item > -90) and (item < 90)

def _validate_ra(item):
    "Ensure a valid right ascention value."
    return (item > 0) and (item < 360)


class CreateConfig(
        PipelineTaskConfig,
        pipelineConnections=CreateConnections):
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




class CreateTask(PipelineTask):
    """Prepare inputs for scene modeling photometry of SNIa candidates.
    """

    _DefaultName = "createSceneModelingInputs"
    ConfigClass = CreateConfig

    def __init__(self, **kwargs):
        super().__init__(**kwargs)   # initialise the mother class PipelineTask
        self.log = logging.getLogger(__name__)   # initialise a Python logger for emitting follow-up messages (info, warning, debug)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, visitImage, visitSummary):
        """
        """
        
        """
        repo = "dp2_prep"
        collection_stage4 = "LSSTCam/runs/DRP/DP2/v30_0_6/DM-53881/stage4"
        butler_stage4 =  Butler(repo, collections=collection_stage4)

        try:
            visitImage = butler_stage4.get('visit_image', visit=visit, detector=det)
        except Exception as e:
            "blablabla"
        """

        wcs = visitImage.getWcs()
        psf = visitImage.getPsf()
        photoCalib = visitImage.getPhotoCalib()
        bbox = visitImage.getBBox()
        visitInfo = visitImage.getInfo().getVisitInfo()
        mjd = visitInfo.getDate().get()  # MJD as float

        # Determine ON/OFF for this visit
        # ON = SN is present
        # OFF = no SN expected
        on_flag = self._compute_on_off(mjd)

        cutouts = []
        psfs = []
        wcss = []
        phot_ratios = []
        positions = []
        on_off = []

        for ra, dec in zip(self.config.ra, self.config.dec):
            sky_coord = geom.SpherePoint(ra, dec, geom.degrees)  # convert SN candidate coords in lsst.geom object

            # Check that the SN falls within this detector
            try:
                pixel_coord = wcs.skyToPixel(sky_coord)
            except Exception as e:
                "blablabla"

            pixel_point = geom.Point2I(pixel_coord)

            # 1- Data cutout
            cutout = self._extract_cutout(visitImage, pixel_coord)
            if cutout is None:
                self.log.warning(
                    f"Could not extract cutout for SN at (ra={ra}, dec={dec}).")
                continue

            # 2- PSF at SN position
            psf_image = self._compute_psf(psf, pixel_coord)

            # 3- PSF derivative

            # 4- WCS

            # 5- Photometric ratio (Ri)
            phot_ratio = self._compute_phot_ratio(photoCalib, pixel_coord)
            

            # Append results
            cutouts.append(cutout)
            psfs.append(psf_image)
            wcss.append(sn_wcs)
            phot_ratios.append(phot_ratio)
            positions.append((ra, dec))
            on_off.append(on_flag)

        return Struct(sceneModelingInputs=result_dict)


    # Private

    def _compute_on_off(self, mjd):
        """ Return True if the SN is expected to be present at this MJD.
        A visit is ON if mjd is between the (min, max) pair in the config.
        It is OFF otherwise (before the explosion or after the SN has faded).
        """
        for mjd_min, mjd_max in zip(self.config.zero_flux_min_mjd, self.config.zero_flux_max_mjd):
            if (mjd_min <= mjd) and (mjd <= mjd_max):
                return True
        return False

    def _extract_cutouts(self, visitImage, sky_coord):
        """Extract a square cutout array centred on pixel_coord (SpherePoint object).
        """
        img_data = visitImage.image.array
        bbox = visitImage.getBBox()
        wcs_obj = visitImage.getWcs()

        pixel_point = wcs_obj.skyToPixel(sky_coord)
        x_loc = int(pixel_point.x - bbox.getMinX())
        y_loc = int(pixel_point.y - bbox.getMinY())

        half_size = self.config.cutout_half_size_pixels
        y0, y1 = y_loc - half_size, y_loc + half_size + 1
        x0, x1 = x_loc - half_size, x_loc + half_size + 1

        cutout = img_data[y0:y1, x0:x1]

        # Gestion des bords
        height, width = img_data.shape
        if (x0 < 0) or
        
        return cutout

    def _computePsf(self, exposure, pixelCenter, srcId=None):
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

    def _compute_psf_derivative():
        return

    def _compute_phot_ratio(self, photoCalib, pixel_coord):
        
        return