===========================================
EDGES Calibration File Structure Definition
===========================================
This document is based on
`memo #113 "Receiver calibration procedure document" <http://loco.lab.asu.edu/loco-memos/edges_reports/tom_20180523_Calibration_Steps.pdf>`_,
by Tom Mozden. It *defines* the directory structure and file naming
conventions used for calibration of EDGES receivers.
Calibration observations that do *not* conform to this structure are considered broken,
and should be either fixed or removed.
This standard is enforced by the code `edges-io <https://github.com/edges-collab/edges-io>`_.

Format
------
* Root directory name format
  + ReceiverXX_YYYY_MM_DD_LLL_to_HHH_MHz
  + XX----------Receiver version number (01,02,03)
  + LLL---------Start frequency in MHz
  + HHH---------Stop frequency in MHz
  + YYYY_MM_DD--Calibration start date
      - eg.---------Receiver03_2019_040_to_200_MHz
* Root directory consists of three subdirectories depending on temperature of receiver
    + 15C---------Calibration is done with receiver temperature set to 15 deg Celsius
    + 25C---------Calibration is done with receiver temperature set to 25 deg Celsius
    + 35C---------Calibration is done with receiver temperature set to 35 deg Celsius
* Each Sub directory consists of three sub folders
    + Resistance--Consists of temperature value from thermistor in .csv format
    + S11---------Consists of VNA measurement in sub folders for different loads and receiver.
     - ReceiverReadingXX---------------XX represents the number of runs with same settings
            ReceiverReadingYY.s1p---S11 of receiver, YY is the number of runs
            ShortYY.s1p-------------S11 of internal short, YY is the number of runs
            OpenYY.s1p--------------S11 of internal open, YY is the number of runs
            MatchYY.s1p-------------S11 of internal Matched, YY is the number of runs
      - AntSimX-------------------------X represents different antenna simulators 1,2,3,4
            ExternalYY--------------S11 of AntsimX, YY is the number of runs
            ShortYY.s1p-------------S11 of internal short, YY is the number of runs
            OpenYY.s1p--------------S11 of internal open, YY is the number of runs
            MatchYY.s1p-------------S11 of internal Matched, YY is the number of runs
     - SwitchingstateXX----------------XX represents the number of runs with same settings
            ShortYY.s1p-------------S11 of internal short, YY is the number of runs
            MatchYY.s1p-------------S11 of internal Matched, YY is the number of runs
            OpenYY.s1p--------------S11 of internal open, YY is the number of runs
            ExternalMatchYY.s1p-----S11 of External standard Matched, YY is the number of runs
            ExternalShortYY.s1p-----S11 of External standard Short, YY is the number of runs
            ExternalOpenYY.s1p------S11 of External standard Open, YY is the number of runs
     - LongCableShort---------------Long cable is connected to short
            ExternalYY--------------S11 of Long cable with short connected, YY is the number of runs
            ShortYY.s1p-------------S11 of internal short, YY is the number of runs
            OpenYY.s1p--------------S11 of internal open, YY is the number of runs
            MatchYY.s1p-------------S11 of internal Matched, YY is the number of runs
     - LongCableOpen-------------------Long cable is connected to open
            ExternalYY--------------S11 of Long cable with open connected, YY is the number of runs
            ShortYY.s1p-------------S11 of internal short, YY is the number of runs
            OpenYY.s1p--------------S11 of internal open, YY is the number of runs
            MatchYY.s1p-------------S11 of internal Matched, YY is the number of runs
  - Ambient-------------------------
            ExternalYY--------------S11 of Long cable with open connected, YY is the number of runs
            ShortYY.s1p-------------S11 of internal short, YY is the number of runs
            OpenYY.s1p--------------S11 of internal open, YY is the number of runs
            MatchYY.s1p-------------S11 of internal Matched, YY is the number of runs
     - HotLoad-------------------------
            ExternalYY--------------S11 of Long cable with open connected, YY is the number of runs
            ShortYY.s1p-------------S11 of internal short, YY is the number of runs
            OpenYY.s1p--------------S11 of internal open, YY is the number of runs
            MatchYY.s1p-------------S11 of internal Matched, YY is the number of runs
  + Spectra----- Consists of digitizer output for all loads, there could be multiple files for each load.
      - file name will be in the format Loadname_YYYY_DDD_HH.acq where YYYY is the year DDD is the day and HH is the Hour.
   | Load Name | Original name format | Renaming format |Notes |
         | --- | --- | --- | ---|
         |  Hotload |  2017_139_19.acq |  HotLoad_2017_139_19.acq | |
         | Ambient  | 2017_141_17.acq  |  Ambient_2017_141_17.acq | |
         | AntSimX  |  2017_143_12.acq |  AntsimX_2017_143_12.acq | X is the antenna simulator version |
         | LongCableShort  |  2017_146_00.acq | LongCableShort_2017_146_00.acq  | |
         |  LongCableOpen | 2017_148_00.acq  |  LongCableOpen_2017_148_00.acq | |
