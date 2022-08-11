# read .stm file, return data structure of the first line. 
# Returns the data dictionary with moments, impar, noisepar, starpar
# see cube2.py for the data structure:
# moments = {'var': power.tolist(), 'cov': cov.tolist(),'rnoise': drnoise[0],'rvar': meanrvar,'mcoef':mcoef.tolist()}
# data = {'image': {'impar': impar, 'noisepar': noisepar}, 'moments': moment
# cubepar = {'nx': nx,'nz':nz,'texp':texp,'gain':gain,'date':DATE,'star':star}

import sys
from os.path import splitext
import numpy as np
from scipy import optimize
from ringss_restore_profile.restore_profile.star import Star
from ringss_restore_profile.restore_profile.core import read_json


class ReadMoments(object):
    def __init__(self, ringss_parameters_filename: str):
        self.parameters = read_json(ringss_parameters_filename)
        self.star_catalog = read_json(self.parameters["profrest"]["starcat"])
        self.weights = read_json(self.parameters["profrest"]["weightsfile"])
        self.zmat = read_json(self.parameters["profrest"]["zmat"])
        self.star = Star(ringss_parameters_filename)

    def __format_profile(self, date: str, hr: str, zen: float, flux: int, profile: dict) -> str:

        s = date + ',' + hr + ',{:.2f}'.format(zen)
        s = s + ',' + str(flux)
        s = s + ',  {:.3f}'.format(profile['see2']) + ',{:.3f}'.format(profile['see']) + ',{:.3f}'.format(
            profile['fsee'])
        s = s + ',{:.2f}'.format(profile['wind']) + ',{:.3f}'.format(profile['tau0']) + ',{:.3f}'.format(
            profile['theta0'])
        s = s + ',  {:.3f}'.format(profile['totvar']) + ',{:.3f}'.format(profile['erms'])
        for x in profile['prof']:
            s += ",{:.2f}".format(x)
        return s

    def __parse_line(self, moments_line: str) -> dict:
        moments_list = moments_line.split(",")

        cube_par = {'nx': int(moments_list[6]), 'nz': int(moments_list[5]), 'texp': float(moments_list[2]) * 1e-3,
                    'gain': float(moments_list[4]), 'date': moments_list[0], 'star': moments_list[1].strip()}

        mcoef = 20  # fixed for now

        k0 = 9
        impar = {'backgr': moments_list[k0], 'flux': moments_list[k0 + 1], 'fluxvar': moments_list[k0 + 2],
                 'meanrad': moments_list[k0 + 3], 'rwidth': moments_list[k0 + 4],
                 'xc': moments_list[k0 + 5], 'yc': moments_list[k0 + 6],
                 'xcvar': moments_list[k0 + 7], 'ycvar': moments_list[k0 + 8],
                 'coma': moments_list[k0 + 9], 'angle': moments_list[k0 + 10], 'contrast': moments_list[k0 + 11]}
        noisepar = moments_list[k0 + 12:k0 + 16]

        m = mcoef + 1
        k1 = k0 + 16
        var = moments_list[k1:k1 + m]
        cov = moments_list[k1 + m:k1 + 2 * m]
        k2 = k1 + 2 * m
        rnoise = moments_list[k2]
        rrms = moments_list[k2 + 1]
        meancoef = moments_list[k2 + 2:k2 + 16]
        momentos = {'var': var, 'cov': cov, 'rnoise': rnoise, 'rvar': rrms, 'mcoef': meancoef}
        data_dict = {'image': {'impar': impar, 'noisepar': noisepar}, 'moments': momentos, 'cube': cube_par}

        return data_dict

    # read & process file with averaged moments only
    # based in read_stm()
    # data_filename, starcat, weights and zmatrix are inputs
    def read_from_file(self, data_filename: str):
        try:
            fh_in = open(data_filename, 'r')  # input .stm file
        except FileNotFoundError as err:
            print(err)
            return None

        filename, file_extension = splitext(data_filename)
        try:
            fh_out = open(filename + '.prof', 'w')  # output profile file
        except FileNotFoundError as err:
            print(err)
            return None

        for line in fh_in:
            data = self.__parse_line(line)

            # star_par = get_starpar(date, hr)
            starpar = self.star.get_starpar(data["cube"]["date"], data["cube"]["star"])

            data["star"] = starpar

            profile = self.restore(data)  # returns profile or None

            if profile:
                print(data['cube']['date'])
                s = self.__format_profile(data['cube']['date'], data['cube']['star'], data['star']['zen'],
                                          data["image"]["impar"]["flux"], profile)
                fh_out.write(s + '\n')

        fh_out.close()
        fh_in.close()
        print('Processing finished!')

    # it is a short version of read_from_file()
    # intended to use directly by RINGSS
    # restore profile from moments & cube data (camera)
    # moments is actually the average moments of all cubes
    # TODO:
    #       I need to figure out how to average dictionaries
    def restore_from_line(self, moments, cube_params):
        data = moments

        # First add cube parameters to data dictionary
        data["cube"] = cube_params

        # get star parameters
        starpar = self.star.get_starpar(cube_params["cube"]["date"], cube_params["cube"]["star"])
        # Add to data dictionary
        data["star"] = starpar

        profile = self.restore(data)  # returns profile or None

        if profile:
            s = self.__format_profile(data['cube']['date'], data['cube']['star'], data['star']['zen'],
                                      data["image"]["impar"]["flux"], profile)
            print(data['cube']['date'])
            print(s)
            return profile

    # Profile restoration. Input: dictionaries of data
    # Output: dictionary of profile parameters
    # It comes from restore() in profrest5.py
    # data is a dictionary. Its structure is the same as calc_moments output plus star parameters and cube parameters
    def restore(self, data: dict) -> dict:
        #    print(data["image"]["impar"])
        var = data["moments"]["var"]  # variance of a-coefficients
        cov = data["moments"]["cov"]  # covariance of a-coefficients
        # impar = data["image"]["impar"]  # ring parameters
        # noisepar = data["image"]["noisepar"]  # noise parameters
        starpar = data["star"]  # star parameters
        zen = starpar["zen"]
        bv = starpar["BV"]
        z0 = self.parameters["profrest"]["zgrid"]  # grid of heights representing the turbulence profile

        z = np.array(self.weights["z"])  # nominal distance grid of weights
        nz = len(self.weights["z"])  # number of layers in weights
        nm = len(self.weights["wt0"][0])  # number of coefficients per layer
        wt = np.reshape(np.array(self.weights["wt0"]), (nz, nm))  # wt.shape = (16,21)
        wt += bv * np.reshape(np.array(self.weights["wtslope"]), (nz, nm))  # must be >0!

        # interpolate weight to the Z0 grid, trim to mmax
        cosz = np.cos(zen / 180 * np.pi)  # cos(z)
        mmax = 15  # hard-coded number of terms to use
        nz0 = len(z0)
        z00 = np.array(z0) / cosz  # stretched array of nominal heights
        wt1 = np.zeros((nz0, mmax))  # interpolated and trimmed weights
        for m in range(0, mmax):
            wt1[:, m] = np.interp(z00, z, wt[:, 1 + m])
        wtsect = np.interp(z00, z, wt[:, 0])  # sector-motion weight on z0 grid

        # Check that moments and weights have the same number of coefficients
        mcoef = len(var)  # 21 for mmax=20
        if mcoef != nm:
            print("Error! Mismatching number of coefficients: ", mcoef, nm)
            return

        # noise bias, see allcubes5.pro => noisecubes
        gain = data["cube"]["gain"]  # electrons per ADU
        eladu = 3.60 * pow(10, -gain / 200)
        noisepar = data["image"]["noisepar"]  # list of 4 numbers
        fluxadu = float(data["image"]["impar"]["flux"])
        #    print(fluxadu)
        flux = eladu * fluxadu  # flux in electrons

        anoise = float(noisepar[0]) / flux + float(noisepar[1]) * pow(self.parameters["telescope"]["ron"] / flux,
                                                                      2)  # noise variance of a-coef
        rnoise = 2. * (float(noisepar[2]) / flux + float(noisepar[3]) * pow(self.parameters["telescope"]["ron"] / flux,
                                                                            2) / 8)  # noise on radius, pix^2

        #    anoise = noisepar[0]/flux + noisepar[1]*pow(par["telescope"]["ron"]/flux,2) # noise variance of a-coef
        #    rnoise = 2.*(noisepar[2]/flux + noisepar[3]*pow(par["telescope"]["ron"]/flux,2)/8) # noise on radius, pix^2

        # Select max frequency used in the restoration
        var1 = np.array(var[1:mmax + 1], float) - anoise
        var1 *= (var1 > 0)  # prevent negative
        cov1 = np.array(cov[1:mmax + 1], float)
        rho = cov1 / var1
        varcorr = var1 / (0.8 + 0.2 * rho)  # correction for finite exposure time
        totvar = np.sum(np.array(var, float))  # full scintillation power, incl. noise

        # Z-matrix correction. Avoid negative values!
        z1 = np.array(self.zmat)
        ncol = z1.shape[1]
        var2 = np.array(var[1:ncol + 1], float)  # indices used for correction
        varz = varcorr / (1. + np.dot(z1[0:mmax, :], var2))  # correct for saturation

        # weighted nnls. Weight is proportional to 1/var (before noise subtraction, non-negative)
        varwt = np.power(np.array(var[1:mmax + 1], float), -1)
        # varwt = np.power(np.array(var[1:mmax+1],float), -0.5)
        a2 = np.transpose(wt1.copy())  # matrix of (mmax-1,nz) dimension
        for i in range(0, mmax):  # weighted system matrix
            a2[i, :] *= varwt[i]
        varz2 = varz * varwt  # weighted right-hand vector
        prof, resvar = optimize.nnls(a2, varz2)  # turbulence integrals in  [m^1/3] in z0 layers

        varmod = np.dot(a2, prof) / varwt
        erms = np.std(1. - varmod / varz)
        print("RMS residual: {:.3f}".format(erms))

        prof1 = prof * cosz  # zenith-corrected integrals
        jtot = np.sum(prof1)  # turbulence integral
        seeconst = 6.83e-13  # integral for 1" seeing at 500nm
        see = pow(jtot / seeconst, 0.6)  # seeing in arcseconds
        jfree = np.sum(prof1[2:nz0 - 1])  # start at 0.5km
        fsee = pow(jfree / seeconst, 0.6)  # FA seeing in arcseconds
        prof1 = prof * cosz * 1e13  # scaled

        # Compute seeing from radius var.
        d = self.parameters["telescope"]["D"]
        pixel = self.parameters["telescope"]["pixel"]
        lameff = self.weights["lameff"]
        lam0 = lameff[0] + bv * (lameff[1] - lameff[0])  # effective wavelength for the star color
        rvar = float(data["moments"]["rvar"]) - rnoise  # radius variance in pix^2, noise-subtracted
        lamd = lam0 / d * 206265.  # lambda/D in arcseconds
        rvarnorm = rvar * pow(pixel / lamd, 2)  # variance in (lam/D^2) units

        wcoef = np.sum(wtsect * prof) / np.sum(prof)  # profile-adjusted weight of sector variance
        jtot2 = rvarnorm / wcoef / 4  # turbulence integral, m^1/3. Explain factor 4!
        see2 = pow(jtot2 * cosz / seeconst, 0.6)  # seeing at zenith, arcsecs
        see2 *= 1. / (1. - 0.4 * totvar)  # saturation correction
        print("Seeing (sect,tot,FA): {:.3f} {:.3f} {:.3f}".format(see2, see, fsee))

        # Wind measurement
        texp = 1e-3  # hard-coded exposure time
        delta = (var1 - cov1) * (texp ** (-2))
        delta = delta * (delta > 0)
        ucoef = np.array(self.weights["ucoef0"]) + bv * np.array(self.weights["ucoefslope"])
        umm = np.array(self.weights["umm"], int) - 1
        delta = np.array(delta[umm])
        # print(delta)  # debugging
        v2mom = np.sum(ucoef * delta)
        jwind = np.sum(prof[2:nz0])  # exclude ground layer, the speed is not zen-corrected
        #    jwind = np.sum(prof[1:nz0])  # exclude ground layer, the speed is not zen-corrected
        v2 = pow(v2mom / jwind, 0.5)  # use profile uncorrected for zenith
        print("Wind speed [m/s]: {:.3f}".format(v2))

        # tau_0 at 500nm
        r0 = pow(6.680e13 * jwind, -0.6)
        tau0 = 310. * r0 / v2  # in [ms]
        # isoplanatic angle at 500nm at zenith
        r0 = 0.101 / see
        tmp = pow(np.array(z0), 1.6667)
        heff = pow(np.sum(tmp * prof1) / np.sum(prof1), 0.6)
        theta0 = 205265. * 0.31 * r0 / heff
        # Output dictionary
        profile_dict = {"z0": z0, "prof": prof1.tolist(), "see": see, "fsee": fsee, "see2": see2, "wind": v2,
                        "erms": erms,
                        "totvar": totvar, "tau0": tau0, "theta0": theta0}

        return profile_dict

    # Calculate moments from given data cube
    # It comes from moments() in cube2.py
    def calc_moments(self, data_cube):
        # print("inside Moments()")
        nz, ny, nx = data_cube.shape  # 2000 64 64

        nstart = 50  # initial average of first 50 frames
        imav = np.average(data_cube[0:nstart], axis=0)  # average nstart frames

        i = np.indices((nx, nx))  # to create xx and yy vectors
        yy = i[0] - nx / 2
        xx = i[1] - nx / 2

        r = np.sqrt(np.square(xx) + np.square(yy))  # distance from center: plt.imshow(xx, cmap='Greys')
        phi = np.arctan2(yy, xx)  # 2D array of phase # plt.imshow(phi, cmap='Greys')

        # ## Prepare initial wide masks
        m = 20  # max order of angular signals
        nsect = 8  # number of sectors for radius calculation
        rwt = np.zeros((nsect, nx, nx))  # for radius calculation
        fwt = np.zeros((nsect, nx, nx))  # fluxes in the sectors
        sect = 2. * np.pi / nsect  # sector width in radians

        for j in range(0, nsect):
            sector = (phi >= (sect * (j - nsect / 2))) & (phi < (sect * (j + 1 - nsect / 2)))
            fwt[j] = sector
            rwt[j] = sector * r

        test = np.zeros((nx, nx))  # test correctness of radial masks
        for j in range(0, nsect):
            test += rwt[j] * (j + 1)

        phisect = sect * (np.arange(nsect, dtype=float) - nsect / 2 + 0.5)  # sector angles
        xsect, ysect = np.cos(phisect), np.sin(phisect)

        backgr = np.median([imav[:, 0], imav[:, nx - 1]])  # left and right columns,  scalar
        tmp = imav - backgr

        itot = np.sum(tmp)  # total flux and centroids
        xc = np.sum(tmp * xx) / itot
        yc = np.sum(tmp * yy) / itot

        imavcent = np.roll(tmp, (int(-xc), int(-yc)), (1, 0))  # crude shift by integer pixel number

        radii = np.zeros(nsect)  # preliminary radii of the sectors, in pixels
        for j in range(0, nsect):
            radii[j] = np.sum(imavcent * rwt[j]) / np.sum(imavcent * fwt[j])

        dx1 = np.sum(radii * xsect) / nsect * 2.3  # accurate x-shift
        dy1 = np.sum(radii * ysect) / nsect * 2.3  # accurate y-shift

        xc += dx1  # more accurate ring center
        yc += dy1

        # FFT sub-pixel shift by -dx1,-dy1
        arg = 2 * np.pi * (dx1 * xx + dy1 * yy) / nx
        imavcent = np.fft.ifft2(np.fft.fft2(imavcent) * np.fft.fftshift(np.cos(arg) + np.sin(arg) * 1j)).real

        # re-compute radii to check the centering, should be similar
        for j in range(0, nsect):
            radii[j] = np.sum(imavcent * rwt[j]) / np.sum(imavcent * fwt[j])
        radpix = np.sum(radii) / nsect

        # Threshold the image at 0.1*max to find the ring width
        tmp = (imavcent - 0.1 * np.max(imavcent))
        tmp *= (tmp > 0)  # plt.imshow(tmp, cmap='Greys')
        radvar = np.sum(tmp * (r - radpix) ** 2) / np.sum(tmp)
        rwidth = pow(radvar, 0.5) * 2.35
        # print("Ring width [pix]: ",rwidth)

        backgr += np.median(imavcent * (r > 1.5 * radpix))  # outside-ring pixels for refined background estimate

        drhopix = 1.5 * rwidth  # mask width, replace 1.5 with parameter value in the future
        ringmask = (r >= radpix - drhopix) * (r <= radpix + drhopix)

        # ### Now build the final matrix of masks
        ncoef = 2 * nsect + 2 * (m + 1)
        maskmat = np.zeros((ncoef, nx * nx))

        for j in range(0, nsect):  # radial masks
            maskmat[j, :] = np.ndarray.flatten(rwt[j] * ringmask)  # image pixels arranged in 1D array
            maskmat[j + nsect, :] = np.ndarray.flatten(fwt[j] * ringmask)

        for j in range(0, m + 1):  # cosine ans sine masks
            tmp = np.cos(phi * j) * ringmask
            if j > 0:
                cwt = tmp - np.sum(tmp) / nx / nx  # remove residual piston
            else:
                cwt = tmp
            tmp = np.sin(phi * j) * ringmask
            swt = tmp - np.sum(tmp) / nx / nx  # remove piston
            maskmat[2 * nsect + j, :] = np.ndarray.flatten(cwt)
            maskmat[2 * nsect + m + 1 + j, :] = np.ndarray.flatten(swt)

        # ### Main loop over the cube
        coef = np.zeros((nz, ncoef))  # prepare the arrays for cube processing
        xcent = np.zeros(nz)  # x-center in each frame [pix]
        ycent = np.zeros(nz)  # y-center [pix]
        rad = np.zeros(nz)
        imav = np.zeros((nx, nx))  # average image
        x0 = xc  # current ring center
        y0 = yc

        for i in range(0, nz):  # process full cube
            tmp = data_cube[i] - backgr
            arg = 2 * np.pi * (x0 * xx + y0 * yy) / nx  # FFT centering
            tmp = np.fft.ifft2(np.fft.fft2(tmp) * np.fft.fftshift(np.cos(arg) + np.sin(arg) * 1j)).real
            imav = imav + tmp
            c = np.dot(maskmat, np.ndarray.flatten(tmp))
            c[0:nsect] = c[0:nsect] / c[nsect:2 * nsect]
            coef[i, :] = c
            radii = c[0:nsect]
            dr = np.sum(radii) / nsect
            dx = np.sum(radii * xsect) / nsect * 2.3
            dy = np.sum(radii * ysect) / nsect * 2.3
            x0 += dx
            y0 += dy
            xcent[i] = x0
            ycent[i] = y0
            rad[i] = dr
            # end of main loop

        # ### Normalization and average parameters of the cube
        imav = imav / nz  # average ring image, save as FITS file
        # hdr = fits.Header()
        # fits.writeto('avimage.fits', imav, hdr) # needs overwrite flag

        flux = np.mean(coef[:, 2 * nsect])  # average flux in the ring [ADU]
        coef[:, 2 * nsect:ncoef] *= 1. / flux  # normalize by the mean flux
        fluxvar = np.std(coef[:, 2 * nsect:ncoef])

        xc = np.mean(xcent)  # mean ring position and its variance
        xcvar = np.std(xcent)
        yc = np.mean(ycent)
        ycvar = np.std(ycent)

        cm = np.mean(coef[:, 2 * nsect + 1])  # mean cosine and sine terms to evaluate coma
        sm = np.mean(coef[:, 2 * nsect + m + 2])
        coma = pow((cm ** 2 + sm ** 2), 0.5)
        angle = 180 / np.pi * np.arctan2(sm, cm)

        contrast = np.zeros(nsect)  # Analog of Strehl ratio in each sector
        for j in range(0, nsect):
            tmp = imav * fwt[j]
            contrast[j] = np.max(tmp) / np.sum(tmp)

        contrast = np.mean(contrast)
        meanrad = np.mean(rad)
        impar = [backgr, flux, fluxvar, meanrad, rwidth, xc, yc, xcvar, ycvar, coma, angle, contrast]  # list object

        tmp = imav / np.sum(imav)  # noise coef. of angular coefficients
        t0 = np.sum(tmp * ringmask)
        noise1 = np.sum(ringmask ** 2 * tmp) / t0
        noise2 = np.sum(ringmask ** 2) / t0

        tmp2 = ringmask * (r - radpix)  # noise coef. of radii
        noise1r = np.sum(tmp * tmp2 ** 2)
        noise2r = np.sum(tmp2 ** 2)
        noisepar = [noise1, noise2, noise1r, noise2r]

        # ###  Calculation of statistical moments
        # Compute differential radius variance
        dr = coef[:, 0:int(nsect / 2)] + coef[:, int(nsect / 2):nsect]  # 4 DIMM-like signals (2000,4)
        drvar = np.var(dr, axis=0)  # 4 variances in pix^2
        meanrvar = np.mean(drvar)

        # radius noise from difference of successive dr values in 4 pairs of opposite sectors
        ddr = dr - np.roll(dr, 1, 0)  # difference (2000,4)
        drnoise = np.var(ddr, axis=0)  # 4 noise variances, pix^2

        # Variance and covariance of angular coefficients
        acoef = coef[:, 2 * nsect:ncoef]  # (2000,42)
        varcoef = np.var(acoef, axis=0)  # 42 variances for m=20
        meancoef = np.mean(acoef, axis=0)
        tmp = acoef * np.roll(acoef, 1, 0)  # shift-1 product
        tmp = tmp[1:nz, :]  # discard first element
        covar = np.sum(tmp, axis=0) / (nz - 1) - meancoef ** 2

        # add cosine and sine variances and covariances
        power = varcoef[0:m + 1] + varcoef[m + 1:2 * m + 2]
        cov = covar[0:m + 1] + covar[m + 1:2 * m + 2]

        # Mean coefficients for aberrations
        mcoef = np.zeros(nsect + 6)  # mean radii and m=1,2,3 cos/sine terms
        mcoef[0:nsect] = np.mean(coef[:, 0:nsect], axis=0)
        mcoef[nsect:nsect + 3] = meancoef[1:4]
        mcoef[nsect + 3:nsect + 6] = meancoef[m + 2:m + 5]

        # format output
        moments = {'var': power.tolist(), 'cov': cov.tolist(), 'rnoise': drnoise[0], 'rvar': meanrvar,
                   'mcoef': mcoef.tolist()}
        data = {'image': {'impar': impar, 'noisepar': noisepar}, 'moments': moments}

        return data, imav


# ### Main module. usage: > python <par> <data>
# ---------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python read_moments.py <par-file.json> <file.stm>")
        sys.exit()

    parameters_filename = sys.argv[1]
    moment = ReadMoments(parameters_filename)

    stm_filename = sys.argv[2]
    print('STM file: ' + stm_filename)

    # Ingest all information needed

    moment.read_from_file(stm_filename)
