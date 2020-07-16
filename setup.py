from cx_Freeze import setup, Executable

executables = [Executable('hrtfmixer.py', base='Win32GUI')]
build_exe_options = {'packages': ['pysofaconventions', 'scipy.spatial', 'matplotlib.pyplot','mpl_toolkits.mplot3d','scipy.signal','numpy','pyaudio','wave','time','pygame'],
'include_files': ['resources/THK_FFHRIR/HRIR_L2354.sofa']}

setup(name='hrtf_mixer',
    version='0.0',
    options={'build_exe': build_exe_options},
    description='Spatialize audio in real time',
    executables=executables)