"""This file contains the PSRFITS template fields. Those fields are defined in
the page :
https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html
"""

main_header = [('HDRVER', '6.1', 'Header version'),
               ('FITSTYPE', 'PSRFITS', 'FITS definition for pulsar data'
                'files'),
               ('DATE', ' ', 'File creation UTC date (YYYY-MM-DDThh:mm:ss)'),
               ('OBSERVER', ' ', 'Observer name(s)'),
               ('PROJID', ' ', 'Project name'),
               ('TELESCOP', ' ', 'Telescope name'),
               ('ANT_X', 0,  '[m] Antenna ITRF X-coordinate (D)'),
               ('ANT_Y', 0,  '[m] Antenna ITRF Y-coordinate (D)'),
               ('ANT_Z', 0,  '[m] Antenna ITRF Z-coordinate (D)'),
               ('FRONTEND', ' ', 'Receiver ID'),
               ('IBEAM', ' ', 'Beam ID for multibeam systems'),
               ('NRCVR', 0, 'Number of receiver polarisation channels'),
               ('FD_POLN', ' ', 'LIN or CIRC'),
               ('FD_HAND', 0, '+/- 1. +1 is LIN:A=X,B=Y, CIRC:A=L,B=R (I)'),
               ('FD_SANG', 0, '[deg] FA of E vect for equal sig in A&B (E)'),
               ('FD_XYPH', 0, '[deg] Phase of A* B for injected cal (E)'),
               ('BACKEND', ' ', 'Backend ID'),
               ('BECONFIG', ' ', 'Backend configuration file name'),
               ('BE_PHASE', 0, '0/+1/-1 BE cross-phase:0 unknown,+/-1 std/rev'),
               ('BE_DCC' , 0, '0/1 BE downconversion conjugation corrected'),
               ('BE_DELAY', 0, '[s] Backend propn delay from digitiser input'),
               ('TCYCLE' , 0, '[s] On-line cycle time (D)'),
               ('OBS_MODE', ' ', '(PSR, CAL, SEARCH)'),
               ('DATE-OBS', ' ', 'UTC date of observation (YYYY-MM-DDThh:mm:ss)'),
               ('OBSFREQ', 0, '[MHz] Centre frequency for observation'),
               ('OBSBW', 0, '[MHz] Bandwidth for observation'),
               ('OBSNCHAN', 0, 'Number of frequency channels (original)'),
               ('CHAN_DM', 0, '[cm-3 pc] DM used for on-line dedispersion'),
               ('PNT_ID', ' ', 'Name or ID for pointing ctr (multibeam feeds)'),
               ('SRC_NAME', ' ', 'Source or scan ID'),
               ('COORD_MD', ' ', 'Coordinate mode (J2000, GALACTIC, ECLIPTIC)'),
               ('EQUINOX', 0, 'Equinox of coords (e.g. 2000.0)'),
               ('RA', ' ', 'Right ascension (hh:mm:ss.ssss)'),
               ('DEC', ' ', 'Declination (-dd:mm:ss.sss)'),
               ('BMAJ', 0, '[deg] Beam major axis length'),
               ('BMIN', 0, '[deg] Beam minor axis length'),
               ('BPA' , 0, '[deg] Beam position angle'),
               ('STT_CRD1', ' ', 'Start coord 1 (hh:mm:ss.sss or ddd.ddd)'),
               ('STT_CRD2', ' ', 'Start coord 2 (-dd:mm:ss.sss or -dd.ddd)'),
               ('TRK_MODE', ' ', 'Track mode (TRACK, SCANGC, SCANLAT)'),
               ('STP_CRD1', ' ', 'Stop coord 1 (hh:mm:ss.sss or ddd.ddd)'),
               ('STP_CRD2', ' ', 'Stop coord 2 (-dd:mm:ss.sss or -dd.ddd)'),
               ('SCANLEN' ,0, '[s] Requested scan length (E)'),
               ('FD_MODE', ' ', 'Feed track mode - FA, CPA, SPA, TPA'),
               ('FA_REQ', 0, '[deg] Feed/Posn angle requested (E)'),
               ('CAL_MODE', ' ', 'Cal mode (OFF, SYNC, EXT1, EXT2)'),
               ('CAL_FREQ', 0, '[Hz] Cal modulation frequency (E)'),
               ('CAL_DCYC', 0, 'Cal duty cycle (E)'),
               ('CAL_PHS', 0, 'Cal phase (wrt start time) (E)'),
               ('CAL_NPHS', 0, 'Number of states in cal pulse (I)'),
               ('STT_IMJD', 0, 'Start MJD (UTC days) (J - long integer)'),
               ('STT_SMJD', 0, '[s] Start time (sec past UTC 00h) (J)'),
               ('STT_OFFS', 0.0, '[s] Start time offset (D)'),
               ('STT_LST', 0, '[s] Start LST (D)')]

subint_header = {'BITPIX': {'value': 8,
                            'comment': 'N/A',
                            'description': ""},
                 'NAXIS': {'value': 2,
                           'comment': '2-dimensional binary table',
                           'description': ""},
                 'NAXIS1': {'value': 0,
                            'comment': 'width of table in bytes',
                            'description': ""},
                 'NAXIS2': {'value': 0,
                            'comment': 'Number of rows in table (NSUBINT)',
                            'description': ""},
                 'PCOUNT': {'value': 0,
                            'comment': 'size of special data area',
                            'description': ""},
                 'GCOUNT': {'value': 1,
                            'comment': 'one data group (required keyword)',
                            'description': ""},
                 #'TFIELDS', 0, 'Number of fields per row', 'description': ""},
                 'EPOCHS': {'value': 'STT_MJD',
                            'comment': 'Epoch convention (VALID, MIDTIME,'
                                       ' STT_MJD)',
                            'description': "Different instruments use different"
                                           " conventions for folding a time"
                                           " series into a pulse profile and"
                                           " for recording the epoch of a"
                                           " sub-integration. This parameter"
                                           " describes one of the following"
                                           " three conventions; if it is not"
                                           " specified in the PSRFITS file,"
                                           " then EPOCHS=STT_MJD is assumed.\n"
                                           "\t*If EPOCHS=VALID, then\n"
                                           "\t\t1.the rising edge of phase"
                                           " bin zero corresponds to an integer"
                                           " turn of the pulsar as predicted by"
                                           " the approximation to the timing"
                                           " model used during folding (e.g."
                                           " tempo2 predictor); and\n"
                                           "\t\t2.the epoch recorded for"
                                           " each sub-integration (typically"
                                           " near the mid-time of the"
                                           " integration) also corresponds to"
                                           " an integer turn of the pulsar.\n"
                                           "\t*If EPOCHS=MIDTIME, then\n"
                                           "\t\t1.the rising edge of phase bin"
                                           " zero corresponds to an integer"
                                           " turn of the pulsar as predicted by"
                                           " the approximation to the timing"
                                           " model used during folding (e.g."
                                           " tempo2 predictor); and\n"
                                           "\t\t2.the epoch recorded for each"
                                           " sub-integration has no relation to"
                                           " pulse phase.\n"
                                           "\t*If EPOCHS=STT_MJD, then\n"
                                           "\t\t1.the rising edge of phase bin"
                                           " zero corresponds to the phase of"
                                           " the pulsar as predicted by the"
                                           " approximation to the timing model"
                                           " used during folding at the epoch"
                                           " defined by the STT_IMJD, STT_SMJD,"
                                           " and STT_OFFS parameters; this is"
                                           " typically the case when the first"
                                           " sample in the time series is"
                                           " folded into the first phase bin."
                                           " Furthermore,\n\t\t2.the epoch"
                                           " recorded for each sub-integration"
                                           " has no relation to pulse phase."},
                 'INT_TYPE': {'value': ' ',
                              'comment': 'Time axis (TIME, BINPHSPERI, BINLNGASC, etc',
                              'description': ""},
                 'INT_UNIT': {'value': ' ',
                              'comment': 'Unit of time axis (SEC, PHS (0-1), DEG)',
                              'description': ""},
                 'SCALE': {'value': ' ',
                           'comment': 'Intensity units (FluxDen/RefFlux/Jansky)',
                           'description': ""},
                 'POL_TYPE': {'value': ' ',
                              'comment': 'Polarisation identifier (e.g., AABBCRCI, AA+BB)',
                              'description': ""},
                 'NPOL': {'value': 0,
                          'comment': 'Nr of polarisations',
                          'description': ""},
                 'TBIN': {'value': 0,
                          'comment': '[s] Time per bin or sample',
                          'description': ""},
                 'NBIN': {'value': 0,
                          'comment': 'Nr of bins (PSR/CAL mode; else 1)',
                          'description': ""},
                 'NBIN_PRD': {'value': 0,
                              'comment': 'Nr of bins/pulse period (for gated data)',
                              'description': ""},
                 'PHS_OFFS': {'value': 0,
                              'comment': 'Phase offset of bin 0 for gated data',
                              'description': ""},
                 'NBITS': {'value': 0,
                           'comment': 'Nr of bits/datum (SEARCH mode data, else 1)',
                           'description': ""},
                 'ZERO_OFF': {'value': 0,
                              'comment': 'Zero offset for SEARCH-mode data',
                              'description': ""},
                 'SIGNINT': {'value': 0,
                             'comment': '1 for signed ints in SEARCH-mode data, else 0',
                             'description': ""},
                 'NSUBOFFS': {'value': 0,
                              'comment': 'Subint offset (Contiguous SEARCH-mode files)',
                              'description': ""},
                 'NCHAN': {'value': 0,
                           'comment': 'Number of channels/sub-bands in this file',
                           'description': ""},
                 'CHAN_BW': {'value': 0,
                             'comment': '[MHz] Channel/sub-band width',
                             'description': ""},
                 'DM': {'value': 0,
                        'comment': '[cm-3 pc] DM for post-detection dedisperion',
                        'description': ""},
                 'RM': {'value': 0,
                        'comment': '[rad m-2] RM for post-detection deFaraday',
                        'description': ""},
                 'NCHNOFFS': {'value': 0,
                              'comment': 'Channel/sub-band offset for split files',
                              'description': ""},
                 'NSBLK': {'value': 0,
                           'comment': 'Samples/row (SEARCH mode, else 1)',
                           'description': ""},
                 'NSTOT': {'value': 0,
                           'comment': 'Total number of samples (SEARCH mode, else 1)',
                           'description': ""}
                }

subint_columns = {'INDEXVAL': {'format': '1D', 'unit': '',
                   'comment': 'Optionally used if INT_TYPE != TIME',
                   'description': "If INT_TYPE is not 'TIME', this column gives"
                                  " the value of the time-like coordinate at"
                                  " the sub-integration centre, expressed in"
                                  " appropriate units (e.g., degrees for"
                                  " longitude)."},
                  'TSUBINT': {'format': '1D', 'unit': 's',
                   'comment':'Length of subintegration',
                   'description':"Duration of sub-integration (or row for"
                                 " search-mode data)"},
                  'OFFS_SUB': {'format': '1D', 'unit': 's',
                   'comment': "Offset from Start of subint centre ",
                   'description': "Time since the observation start at the"
                                  "centre of each sub-integration (or row)." },
                  'LST_SUB': {'format': '1D', 'unit': 's',
                   'comment': "LST at subint centre",
                   'description': "Approximate local sidereal time at the"
                                  " sub-integration centre"},
                  'RA_SUB': {'format': '1D', 'unit': 'deg',
                   'comment': "RA (J2000) at subint centre",
                   'description': "Pointing J2000 Right Ascension at the time"
                                  " of the sub-integration centre. For scanning"
                                  " observations this may change with time." },
                  'DEC_SUB': {'format': '1D', 'unit': 'deg',
                   'comment': "Dec (J2000) at subint centre ",
                   'description': "Pointing J2000 Declination at the time of"
                                  " the sub-integration centre. For scanning"
                                  " observations this may change with time." },
                  'GLON_SUB': {'format': '1D', 'unit': 'deg',
                   'comment': "[deg] Gal longitude at subint centre",
                   'description': "Pointing Galactic longitude at the time of"
                                  " the sub-integration centre. For scanning"
                                  " observations this may change with time." },
                  'GLAT_SUB': {'format': '1D', 'unit': 'deg',
                   'comment': "[deg] Gal latitude at subint centre",
                   'description':  "Pointing Galactic latitude at the time of"
                                   " the sub-integration centre. For scanning"
                                   " observations this may change with time." },
                  'FD_ANG': {'format': '1E', 'unit': 'deg',
                   'comment': "[deg] Feed angle at subint centre",
                   'description': "Angle of the feed reference (normally the"
                                  " A or X probe) with respect to the telescope"
                                  " zenith meridian." },
                  'POS_ANG': {'format': '1E', 'unit': 'deg',
                   'comment': "[deg] Position angle of feed at subint centre",
                   'description': "Angle of the feed reference (normally the"
                                  " A/X probe) with respect to the celestial"
                                  " meridian. This may be held at a fixed value"
                                  " by adjusting the feed angle according to"
                                  " the variation of parallactic angle during"
                                  " an observation (FD_MODE = 'CPA')." },
                  'PAR_ANG': {'format': '1E', 'unit': 'deg',
                   'comment': "[deg] Parallactic angle at subint centre",
                   'description': "Parallactic angle at the time of the"
                                  " sub-integration centre." },
                  'TEL_AZ': {'format': '1E', 'unit': 'deg',
                   'comment': "[deg] Telescope azimuth at subint centre",
                   'description': "Telescope azimuth angle at the time of the"
                                  " sub-integration centre." },
                  'TEL_ZEN': {'format': '1E', 'unit': 'deg',
                   'comment': "[deg] Telescope zenith angle at subint centre",
                   'description': "Telescope zenith angle at the time of the"
                                  "sub-integration centre." },
                  'AUX_DM': {'format': '1D', 'unit': 'CM-3 PC',
                   'comment': "additional DM (ionosphere, corona, etc.)",
                   'description': "An additional time-varying component of"
                                  " dispersion measure to be added to the"
                                  " interstellar component. This allows"
                                  " observed dispersion measures to separated"
                                  " into a constant or slowly varying component"
                                  " and a time-varying (normally solar-system)"
                                  " component." },
                  'AUX_RM': {'format': '1D', 'unit': 'RAD M-2',
                   'comment': "additional RM (ionosphere, corona, etc.)",
                   'description': "An additional time-varying component of"
                                  " rotation measure to be added to the"
                                  " interstellar component. This allows"
                                  " observed rotation measures to separated"
                                  " into a constant or slowly varying component"
                                  " and a time-varying (normally solar-system)"
                                  " component." },
                  'DAT_FREQ': {'format': 'D', 'unit': 'MHz',
                   'comment': "[MHz] Centre frequency for each channel ",
                   'description': "Centre frequency of each channel or"
                                  " sub-band. For data where channels have been"
                                  " summed to form sub-bands, the frequency of"
                                  " each sub-band is the weighted mean of the"
                                  " centre frequencies of the summed channels."},
                  'DAT_WTS': {'format': 'E', 'unit': '',
                   'comment': "Weights for each channel",
                   'description': "Channel weights, in the range 0 - 1."
                                  " Channels may be given zero weight if they"
                                  " are affected by interference or outside"
                                  " the effective bandpass. Normally defaults"
                                  " to 1.0 for original recorded data." },
                  'DAT_OFFS': {'format': 'E', 'unit': '',
                   'comment': "Data offset for each channel",
                   'description': "Value subtracted from data to give zero mean"
                                  " for each channel." },
                  'DAT_SCL': {'format': 'E', 'unit': '',
                   'comment': " Data scale factor (outval=dataval*scl + offs)",
                   'description': "Scale factor used so that the integer table"
                                  " data cover the full available range." },
                  'DATA': {'format': 'I', 'unit': 'Jy',
                   'comment': "Subint data table",
                   'description': "Table containing the digitised data."
                                  " For fold-mode data, 16-bit integers are"
                                  " used to represent the scaled and"
                                  " baseline-subtracted data. Search-mode data"
                                  " are stored as a byte array. If NBIT < 8,"
                                  " samples are packed into bytes with earlier"
                                  " samples in the higher-order bits."}
                  }
