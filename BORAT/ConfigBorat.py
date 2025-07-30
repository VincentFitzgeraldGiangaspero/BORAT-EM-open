from configparser import ConfigParser
import ast
import os
import numpy as np

c = 299792458.0  # Speed of light in m/s


class Read:

    def __init__(self, configFile):

        config = ConfigParser()
        config.read(configFile)

# CFD input
        self.Case = config.get('CFD', 'Case')

        self.CFDsolution = config.get('CFD', 'CFDsolution')
        self.CFDsolutionSurface = config.get('CFD', 'InputFileSurface')

        self.Dimension = config.get('CFD', 'Dimension')
        self.precomputedRefractive = config.getboolean(
            'CFD', 'Precomputed_Refractive')
        self.precomputedGrad = config.getboolean(
            'CFD', 'Precomputed_Gradient')

# Antenna

        self.Frequency = config.getfloat(
            'Antenna', 'Frequency')
        self.Antenna_Aperture = config.getfloat(
            'Antenna', 'Antenna_Aperture')
        self.AntennaAxis = ast.literal_eval(config.get(
            'Antenna', 'AntennaAxis'))
        self.Antenna_Location = ast.literal_eval(config.get(
            'Antenna', 'Antenna_Location'))
        self.Antenna_Direction = config.getint(
            'Antenna', 'Antenna_Direction')
        self.TotalRays = config.getint(
            'Antenna', 'TotalRays')
        self.Rho = config.getfloat(
            'Antenna', 'Rho')
        self.iPol = config.getint(
            'Antenna', 'Polarization')
        self.E0_Amplitude = config.getfloat(
            'Antenna', 'E0_Amplitude')

        self.Lambda = c/self.Frequency

        self.K = 2*np.pi/self.Lambda

# RayTracing

        self.Ionized_solution = config.getboolean(
            'RayTracing', 'Ionized_solution')
        self.Reflections = config.getboolean(
            'RayTracing', 'Reflections')

        # self.TotalRays=self.n_alpha*self.n_theta

# EikonalIntegrator

        self.raysSampling = config.getint(
            'EikonalIntegrator', 'raysSampling')
        self.maxNraySteps = config.getint(
            'EikonalIntegrator', 'maxNraySteps')
        self.maxStep = config.getfloat(
            'EikonalIntegrator', 'maxStep')
        self.firstStep = config.getfloat(
            'EikonalIntegrator', 'firstStep')
        self.CutoffFlag = config.getboolean(
            'EikonalIntegrator', 'Cutoff')
        self.CutoffValue = config.getfloat(
            'EikonalIntegrator', 'CutoffValue')
        self.Tolerance = config.getfloat(
            'EikonalIntegrator', 'Tolerance')
        self.method = config.get(
            'EikonalIntegrator', 'method')

# Aperture Integration

        self.tstart_obs = config.getint(
            'Aperture', 'tstart_obs')
        self.tstop_obs = config.getint(
            'Aperture', 'tstop_obs')
        self.delt = config.getfloat(
            'Aperture', 'delt')

        self.pstart_obs = config.getint(
            'Aperture', 'pstart_obs')
        self.pstop_obs = config.getint(
            'Aperture', 'pstop_obs')
        self.delp = config.getfloat(
            'Aperture', 'delp')


# Parallelization

        self.Parallel = config.getboolean(
            'Parallelization', 'Parallel')

        self.nProcess = config.getint(
            'Parallelization', 'Processors')

        self.Parallel_Method = config.get(
            'Parallelization', 'Parallel_Method')

# Output

        self.OutputFolder = config.get(
            'Output', 'OutputFolder')

        self.raySolutionFileTecplot = config.get(
            'Output', 'raySolutionFileTecplot')

        self.raySolutionFilePLK = config.get(
            'Output', 'raySolutionFilePKL')

        if not os.path.exists(self.OutputFolder):
            os.makedirs(self.OutputFolder)

# TecIO

        self.tecioPath = config.get(
            'TecIO', 'Tecio_Path')


#########################################################################
