from pysofaconventions import * # For getting measurement locations and FIRs from SOFA files
from scipy.spatial import Delaunay # For triangulating the measurement points, provides data for pathing algorithm
import matplotlib.pyplot as plt # at this point, mainly just debugging
from mpl_toolkits.mplot3d import Axes3D # also debugging
import scipy.signal # fftconvolve of recorded data with FIRs
import numpy as np # math
import pyaudio # streaming and real-time audio processing
import wave # reading and writing wav files
import time # debugging
import pygame # user interface for hrtf measurement location input
# from pygame import gfxdraw
from mixerGraphics import Circle
import filename_dragger

##### CALLBACK
# While the stream is playing, it can processes the next bit of audio data in callback
# Think of it as this: it plays a chunk's worth of audio at a time, and 
# that acts as downtime for processing the next chunk
# This is where we apply the HRTF

def callback(in_data, frame_count, time_info, status):
    # uses the global paz and pel variables calculated using the pygame knob positions
    px = np.sin(paz)*np.cos(pel)*pr
    py = np.cos(paz)*np.cos(pel)*pr
    pz = np.sin(pel)*pr

    pp = np.array((px, py, pz))

    ### Adjacency Walk to find Tetra
    # By keeping currentTetraIndex global, it prevents guarantees a short path. Faster than octree b/c we know the previous state.
    global currentTetraIndex
    i = 0
    while True:
        [g1, g2, g3] = (pp-tetraCoords[currentTetraIndex,3])@Tinv[currentTetraIndex]
        g4 = 1-g1-g2-g3
        gs = [g1, g2, g3, g4]
        if all(g >= 0 for g in gs) or i>=20000:#len(tetraCoords):
            break
        currentTetraIndex = tri.neighbors[currentTetraIndex][gs.index(min(gs))]
        i+=1

    # and get the HRTF associated with pp
    origPosIndex = tri.simplices[currentTetraIndex]
    hrtfA = FIRs[origPosIndex[0],:,:]
    hrtfB = FIRs[origPosIndex[1],:,:]
    hrtfC = FIRs[origPosIndex[2],:,:]
    hrtfD = FIRs[origPosIndex[3],:,:]

    hrtf = hrtfA*gs[0]+hrtfB*gs[1]+hrtfC*gs[2]+hrtfD*gs[3]
    # hrtf = hrtfA
    global hrtfList, dataPrepend
    if not np.array_equal(hrtfList[-1][0], hrtf[0]):
        hrtfList.append(hrtf)
        # if len(hrtfList)>2:
        #     hrtf = np.array(hrtfList[-2])*.99 + hrtf*.01

    # Now read the data, apply the convolution, and output
    # global overlapLeft, overlapRight
    data = wf.readframes(frame_count)
    data_int = np.frombuffer(data, dtype=np.int16)
    data_int = np.concatenate((dataPrepend, data_int))
    dataPrepend = data_int[-overlapAmount:]

    # Convolve with the hrtf
    binaural_left = scipy.signal.fftconvolve(data_int,hrtf[0])
    binaural_right = scipy.signal.fftconvolve(data_int,hrtf[1])

    # Overlap-add of convolution tails
    if len(binaural_left)>0:
        binaural_left = binaural_left[overlapAmount:-overlapAmount]
        binaural_left = binaural_left.astype(np.int16)
        # soundPlot[0].extend(binaural_left)

        binaural_right = binaural_right[overlapAmount:-overlapAmount]
        binaural_right = binaural_right.astype(np.int16)
        # soundPlot[1].extend(binaural_right)

    # Interleave into stereo byte array
    binaural = np.empty((binaural_left.size + binaural_right.size,), dtype=np.int16)
    binaural[0::2] = binaural_left
    binaural[1::2] = binaural_right
    data = binaural[:CHUNK*2].tobytes()
    recording.append(data)
    return (data, pyaudio.paContinue)

# Everything in here is run once
##### SETUP
soundPlot = [[],[]] # Get the waveform of the output as an int array
hrtfList = [[0,0]] # Get the hrtf profile
if True:
    ##### HRTF Setup
    def setupHRTF():
        # Paths to sofa files
        # You can download these at https://www.sofaconventions.org/mediawiki/index.php/Files
        # Code should work for any spherically recorded HRTFs
        # This one is included in the resources folder of the repository
        folderPath = 'resources/THK_FFHRIR/'
        fileNames = [
            'HRIR_L2354',
            'HRIR_L2354'
        ]

        # Should be noted that the code assumes the software is of the SimpleFreeFieldHRIR standard
        # Not really sure what happens when you use a file that's not of that format...
        sofaFiles = [SOFAFile(folderPath+fileName+'.sofa','r') for fileName in fileNames]

        # Each source position variable is a list of az, el, r
        sourcePositions = np.concatenate([sofaFile.getVariableValue('SourcePosition')
            for sofaFile in sofaFiles])
    
        # Turn az and el into radians
        sourcePositions[:,:2] *= np.pi/180

        # Increasing cull Amount might get rid of the snipping?
        # I think it's happening because the point moves enough to make it sound discontinuous
        # look to see how the pyo package deals with it? Maybe it is just a matter of speeding up
        # how often the azs and els are sampled

        # Maybe subpixel continuity would help too?

        # This actually has a functional purpose outside of less computation. Delaunay triangulation requires
        # a sparse distribution of points to maintain numerical stability. Otherwise you get high aspect ratio
        # tetrahedrons (tets) that can flatten out or have super high magnitude inverses (used in tet searching).
        cullAmount = 3

        # We double the surface of points and offset it outwards in order
        # to create a sheath of delaunay tetrahedrons with good aspect ratios.
        # The mean free path is actually just the approximate point to point distance
        # when you speckle points on a sphere, but mean free path sounded cooler.
        meanFreePath = 4*max(sourcePositions[:,2])/np.sqrt(len(sourcePositions)/cullAmount)
        sourcePositions[len(sourcePositions)//2:,2] += meanFreePath

        maxR = max(sourcePositions[:,2])-meanFreePath/2
        
        # FIR (Finite Impulse Response) in the form (measurement point index, left IR, right IR)
        FIRs = np.concatenate([sofaFile.getDataIR()
            for sofaFile in sofaFiles])

        # Extract all az, el, r. NOTE: These are in the same order as the sourcePositions list
        az = np.array(sourcePositions[:,0])
        el = np.array(sourcePositions[:,1])
        r = np.array(sourcePositions[:,2]) 

        xs = np.sin(az)*np.cos(el)*r
        ys = np.cos(az)*np.cos(el)*r
        zs = np.sin(el)*r

        # # Spacing Check
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.scatter(xs[:500], ys[:500], zs[:500])
        # ax.scatter(xs[2702:3202], ys[2702:3202], zs[2702:3202])
        # ax.set_xlim([-4,4])
        # ax.set_ylim([-4,4])
        # ax.set_zlim([-4,4])
        # plt.show()
        # aoeu

        # points is now a list of [x,y,z] in the order of sourcePositions
        points = np.array([xs, ys, zs]).transpose()

        
        sourcePositions = sourcePositions[::cullAmount]
        FIRs = FIRs[::cullAmount]
        points = points[::cullAmount]

        # tri is a delaunay object thingy
        tri = Delaunay(points, qhull_options="QJ Pp")

        # Setup barycentric coordinate calculation
        tetraCoords = points[tri.simplices] # List of tetras' coordinates of points

        T = np.transpose(np.array((tetraCoords[:,0]-tetraCoords[:,3],
                    tetraCoords[:,1]-tetraCoords[:,3],
                    tetraCoords[:,2]-tetraCoords[:,3])), (1,0,2))

        def fast_inverse(A):
            identity = np.identity(A.shape[2], dtype=A.dtype)
            Ainv = np.zeros_like(A)
            planarCount=0
            for i in range(A.shape[0]):
                try:
                    Ainv[i] = np.linalg.solve(A[i], identity)
                except np.linalg.LinAlgError:
                    # If there's a flat object, it's going to create an
                    # infinite value for the inverse matrix (det = 0)
                    # Use this value to help debug geometry
                    planarCount += 1
            # print(planarCount)
            return Ainv

        Tinv = fast_inverse(T) # a list of all the barycentric inverses of T, listed in the same order as the tetras in tri and tetraCoords
        return(tetraCoords, Tinv, tri, FIRs, maxR)

    # These are the only parameters from the setup that matter
    tetraCoords, Tinv, tri, FIRs, maxR = setupHRTF()

    # averageFIR = np.zeros((2, FIRs.shape[2]))
    # for FIR in FIRs:
    #     averageFIR += FIR/len(FIRs)

    # fig = plt.figure()
    # ax1 = fig.add_subplot(311)
    # ax2 = fig.add_subplot(312)
    # ax3 = fig.add_subplot(313)
    # ax1.plot(averageFIR[0])
    # ax2.plot(averageFIR[1])
    # ax3.plot(FIRs[0,0])
    # plt.show()


    # If chunk size is too high, then movement will not be smooth because
    # every chunk is played at a particular az/el/r.
    # If chunk size is too low, then your computer will not be able to 
    # process the audio fast enough.
    CHUNK = 30

    # This is the overlap value of the convolution. When you convolve
    # the FIR with the audio signal, it will have a tail that is this long.
    overlapAmount = FIRs.shape[2]-1
    dataPrepend = np.zeros(overlapAmount)

    # Initializing the current tetrahedron index for the callback function
    # By starting the search at our last known index, we greatly reduce the
    # search time required b/c we're probably pretty close to where we were
    # a moment ago.
    currentTetraIndex = 0

    # Roughly half of head diameter
    minR = 0.075

    ##### Recording Setup
    dragger = filename_dragger.dragStarter()
    file = dragger.fileName[0]
    # filePath = 'C:/Users/Alex/Videos/ASMR/'
    # fileName = 'testAudioFile.wav'
    wf = wave.open(file, 'rb')
    p = pyaudio.PyAudio()

    # This is the .wav file output
    recording = []

    sourcePath = []
    # sourcePath = [[np.pi/2,0,2,2],
    #         [-np.pi/2,0,2,2]]

    # #### Recording Viz Setup
    # # A way to see what waveforms are coming up
    # wfviz = wave.open(fileName, 'rb')
    # CHUNKviz = 1024*100
    # dataviz = wfviz.readframes(CHUNKviz)
    # dataviz_int = np.frombuffer(dataviz, dtype=np.int16)
    # plt.plot(dataviz_int[::100])
    # plt.show()

    running = True

    azimuth = 0
    elevation = 0
    dist = maxR
    paz = azimuth
    pel = elevation
    pr = dist

    playing = False

    stream = p.open(
        format=p.get_format_from_width(wf.getsampwidth()),
        channels=2,#wf.getnchannels(),
        rate=wf.getframerate(),
        output=True,
        frames_per_buffer=CHUNK,
        start = False,
        stream_callback=callback)

##### GUI RUNNING ####
def gameGUI(sourcePath=[]):
    global azimuth, elevation, dist, paz, pel, pr, playing

    sourcePath = sourcePath
    posIn = 0
    scrubTime = 0
    limitTime = 0
    startTime = 0
    downTime = 0
    pauseStart = 0

     ##### GUI SETUP
    SCREEN_WIDTH = 1000
    SCREEN_HEIGHT = 1000

    BG_COLOR = (40, 41, 35) # Dark Grey
    color_R = (249, 36, 114) # RED
    color_B = (103, 216, 239)
    color_G = (166, 226, 43)
    color_O = (253, 150, 34)
    color_W = (248, 248, 239)
    color_P = (172, 128, 255)
    color_LG = (116, 112, 93)

    # I leave FPS off b/c if you keep the gui at a particular fps then
    # the audio processing will have to wait until the frame is done
    # to play the next bit of audio. This makes it choppy, and always will
    # be unless you can run at 44,100 fps, and I promise you your TN can't
    # do that, no matter what reddit told you.
    FPS = 60
    clock = pygame.time.Clock()

    pygame.init()
    window = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT), pygame.RESIZABLE)

    def text_objects(text, font, color):
        textSurface = font.render(text, True, color)
        return textSurface, textSurface.get_rect()

    # PLOTS
    orig_x = int(SCREEN_WIDTH/2)
    orig_y = int(SCREEN_HEIGHT*0.4)
    dispRad = int(SCREEN_HEIGHT/4)
    polPlot = Circle(orig_x, orig_y, color_LG, dispRad)

    v_spacing = 50
    azThick = 10
    azPlot_r = SCREEN_HEIGHT//10
    azPlot_x = polPlot.pos[0]
    azPlot_y = polPlot.pos[1]+polPlot.radius+azPlot_r+v_spacing
    azPlot = Circle(polPlot.pos[0], polPlot.pos[1]+polPlot.radius+azPlot_r+v_spacing, color_LG, azPlot_r)

    elPlot = pygame.Rect(int(orig_x*5/3)-2, orig_y-dispRad/2, 4, dispRad)
    rPlot = pygame.Rect(int(orig_x*5/3-dispRad/2)-2, orig_y+dispRad, dispRad, 4)

    # CONTROLS
    polCursor = Circle(orig_x, int(orig_y-90/140*dispRad), color_B, 10)
    azCursor = Circle(orig_x, azPlot_y-azPlot_r, color_B, 10)
    elCursor = Circle(elPlot.centerx, int(elPlot.top+90/140*elPlot.height), color_B, 10)
    rCursor = Circle(int(rPlot.left + (dist-minR)/(maxR-minR)*rPlot.width), 
                        rPlot.centery, color_B, 10)
    cursorList = [polCursor, azCursor, elCursor, rCursor]

    activeCursor = None

    pressing_List = []

    screenChange = False

    def updateGraphics():
        window.fill(BG_COLOR)

        if playing:
            largeText = pygame.font.SysFont('lucida console',20)
            TextSurf, TextRect = text_objects('Press [SPACE] to stop recording...', largeText, color_O)
            TextRect.left = SCREEN_WIDTH/8/2
            TextRect.top = SCREEN_HEIGHT/8/2
            window.blit(TextSurf, TextRect)

            largeText = pygame.font.SysFont('lucida console',40)
            TextSurf, TextRect = text_objects('RECORDING...', largeText, color_R)
            TextRect.left = SCREEN_WIDTH/8/2
            TextRect.top = SCREEN_HEIGHT/8
            window.blit(TextSurf, TextRect)

        if not playing:
            largeText = pygame.font.SysFont('lucida console',20)
            TextSurf, TextRect = text_objects('Press [SPACE] to record...', largeText, color_B)
            TextRect.left = SCREEN_WIDTH/8/2
            TextRect.top = SCREEN_HEIGHT/8/2
            window.blit(TextSurf, TextRect)

        ## Polar Plot
        polPlot.draw_ring(window, 5)

        # Polar Plot Els
        for gradation in [-3, 0, 3, 6]:
            thick = 1
            if gradation == 0:
                thick = 3
            pygame.draw.circle(
                window, color_LG,
                polPlot.pos, int((90-gradation*10)/140*dispRad), thick)

        # Polar Plot Azs
        for gradation in np.linspace(np.pi/2, np.pi*5/2, 9):
            pygame.draw.aaline(
                window, color_W,
                polPlot.pos, (int(polPlot.pos[0]+np.cos(gradation)*polPlot.radius), 
                int(polPlot.pos[1]+np.sin(gradation)*polPlot.radius)))

        polCursor.draw_circle(window)

        ## Az Plot
        azPlot.draw_ring(window, azThick)
        azCursor.draw_circle(window)

        ## El Plot
        pygame.draw.rect(
            window, color_LG,
            elPlot)
        elCursor.draw_circle(window)

        ## R Plot
        pygame.draw.rect(
            window, color_LG,
            rPlot)
        rCursor.draw_circle(window)

        pygame.display.update()

    running = True
    # updateGraphics()
    while running:
        ### EVENTS
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

            elif event.type == pygame.VIDEORESIZE:
                SCREEN_WIDTH = event.w
                SCREEN_HEIGHT = event.h
                window = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT), pygame.RESIZABLE)
                orig_x = int(SCREEN_WIDTH/2)
                orig_y = int(SCREEN_HEIGHT*0.4)
                dispRad = int(SCREEN_HEIGHT/4)
                polPlot = Circle(orig_x, orig_y, color_LG, dispRad)

                v_spacing = 50
                azThick = 10
                azPlot_r = SCREEN_HEIGHT//10
                azPlot_x = polPlot.pos[0]
                azPlot_y = polPlot.pos[1]+polPlot.radius+azPlot_r+v_spacing
                azPlot = Circle(polPlot.pos[0], polPlot.pos[1]+polPlot.radius+azPlot_r+v_spacing, color_LG, azPlot_r)

                elPlot = pygame.Rect(int(orig_x*5/3)-2, orig_y-dispRad/2, 4, dispRad)
                rPlot = pygame.Rect(int(orig_x*5/3-dispRad/2)-2, orig_y+dispRad, dispRad, 4)

                polCursor.pos = (int(np.cos(azimuth-np.pi/2)*polPlot.radius*(90-elevation*180/np.pi)/140)+polPlot.pos[0], 
                                 int(np.sin(azimuth-np.pi/2)*polPlot.radius*(90-elevation*180/np.pi)/140)+polPlot.pos[1])

                screenChange = True
                
            if not sourcePath:
                # Test for MOUSEBUTTONDOWN events
                if event.type == pygame.MOUSEBUTTONDOWN:
                    # User pressed mouse buttons
                    if event.button == 1:
                        pos = pygame.mouse.get_pos()
                        for i, cursor in enumerate(cursorList):
                            if ((pos[0]-cursor.pos[0])**2 + (pos[1]-cursor.pos[1])**2 < cursor.radius**2):
                                # if click on any of the control cursors then 
                                # that's the active cursor now
                                activeCursor = cursor

                                # This keeps the control cursor from snapping
                                # its center to the mouse
                                mouse_x, mouse_y = event.pos
                                offset_x = activeCursor.pos[0] - mouse_x
                                offset_y = activeCursor.pos[1] - mouse_y
                                break

                elif event.type == pygame.MOUSEBUTTONUP:
                    # User released mouse buttons
                    if event.button == 1:
                        activeCursor = None

                elif event.type == pygame.MOUSEMOTION:
                    # Move the plots only, location parameters not changed here
                    # Az, el, and distance are altered in the update section
                    if activeCursor is not None:
                        mouse_x, mouse_y = event.pos
                        activeCursor.pos = (mouse_x + offset_x, mouse_y + offset_y)

                        # Polar Plot
                        if activeCursor == polCursor:
                            pC_rad = polCursor.radDist(polPlot)**0.5
                            if pC_rad > polPlot.radius:
                                pol_xlim = (polCursor.pos[0]-polPlot.pos[0])/pC_rad*polPlot.radius
                                pol_ylim = (polCursor.pos[1]-polPlot.pos[1])/pC_rad*polPlot.radius
                                polCursor.pos = (int(pol_xlim)+polPlot.pos[0],
                                                int(pol_ylim)+polPlot.pos[1])

                        # Az Plot                
                        elif activeCursor == azCursor:
                            aC_rad = azCursor.radDist(azPlot)**0.5
                            if aC_rad != azPlot.radius:
                                az_xlim = (azCursor.pos[0]-azPlot.pos[0])/aC_rad*azPlot.radius
                                az_ylim = (azCursor.pos[1]-azPlot.pos[1])/aC_rad*azPlot.radius
                                azCursor.pos = (int(az_xlim)+azPlot.pos[0],
                                                int(az_ylim)+azPlot.pos[1])

                        # El Plot
                        elif activeCursor == elCursor:
                            elCursor.pos = (elPlot.centerx, elCursor.pos[1])
                            if elCursor.pos[1]>elPlot.bottom:
                                elCursor.pos = (elCursor.pos[0], elPlot.bottom)

                            elif  elCursor.pos[1]<elPlot.top:
                                elCursor.pos = (elCursor.pos[0], elPlot.top)

                        # R Plot
                        elif activeCursor == rCursor:
                            rCursor.pos = (rCursor.pos[0], rPlot.centery)
                            if rCursor.pos[0]>rPlot.right:
                                rCursor.pos = (rPlot.right, rCursor.pos[1])
                            elif  rCursor.pos[0]<rPlot.left:
                                rCursor.pos = (rPlot.left, rCursor.pos[1])

            if event.type == pygame.KEYDOWN:
                key_name = pygame.key.name(event.key)
                pressing_List.append(key_name)
                if key_name == 'space':
                    playing = not playing
                    if playing:
                        if startTime == 0:
                            startTime = time.time()
                        if pauseStart != 0:
                            downTime = downTime + (time.time()-pauseStart)
                        stream.start_stream()
                        if not sourcePath:
                            for cursor in cursorList:
                                cursor.color=color_R
                        else:
                            for cursor in cursorList:
                                cursor.color=color_G
                    else:
                        pauseStart = time.time()
                        stream.stop_stream()
                        for cursor in cursorList:
                            cursor.color=color_B

            elif event.type == pygame.KEYUP:
                key_name = pygame.key.name(event.key)
                pressing_List.remove(key_name)


        ### UPDATES
        for key in pressing_List:
            if key == 'd':
                rCursor.pos = (int(min(rCursor.pos[0]+rPlot.width/100, rPlot.right)), rPlot.centery)
            elif key == 'a':
                rCursor.pos = (int(max(rCursor.pos[0]-rPlot.width/100, rPlot.left)), rPlot.centery)
            # Do implementation of az cursor control at some point
            # Current issue is that the logic needs to be figured out when polCursor is active
            # elif key == 'w':
            #     elCursor.pos = (elPlot.centerx, max(elCursor.pos[1]-1, elPlot.top))
            # elif key == 's':
            #     elCursor.pos = (elPlot.centerx, min(elCursor.pos[1]+1, elPlot.bottom))
            # elif key == 'a':
            #     az_x = azCursor.pos[0]-azPlot.pos[0]
            #     az_y = azCursor.pos[1]-azPlot.pos[1]
            #     d_az = np.array([-az_y/azPlot.radius, az_x/azPlot.radius])*2
            #     azCursor.pos = (int(azCursor.pos[0]-d_az[0]),
            #                     int(azCursor.pos[1]-d_az[1]))
            #     aC_rad = azCursor.radDist(azPlot)**0.5
            #     az_xlim = (azCursor.pos[0]-azPlot.pos[0])/aC_rad*azPlot.radius
            #     az_ylim = (azCursor.pos[1]-azPlot.pos[1])/aC_rad*azPlot.radius
            #     azCursor.pos = (int(az_xlim)+azPlot.pos[0],
            #                     int(az_ylim)+azPlot.pos[1])
            #     print(d_az)
            # elif key == 'd':
            #     azCursor.pos

        # If moved radio, move bars, as all parameters are determined from bars
        if activeCursor == polCursor or screenChange:
            cartesian_x = polCursor.pos[0] - polPlot.pos[0]
            cartesian_y = polCursor.pos[1] - polPlot.pos[1]

            polar_az = np.pi/2-np.arctan2(-cartesian_y, cartesian_x)
            polar_el = 90-np.sqrt(cartesian_x**2 + cartesian_y**2)/dispRad*140

            azCursor.pos = (int(np.cos(polar_az-np.pi/2)*azPlot.radius)+azPlot.pos[0], 
                            int(np.sin(polar_az-np.pi/2)*azPlot.radius)+azPlot.pos[1])
            elCursor.pos = (elPlot.centerx,
                            int(elPlot.bottom-(polar_el+50)/140*elPlot.height))

        if screenChange:
            rCursor.pos = (int((dist-minR)/(maxR-minR)*rPlot.width + rPlot.left),
                            rPlot.centery)
            screenChange = False

        # Set HRTF parameters
        # Everything enters and leaves in radians
        if not sourcePath:
            azimuth = np.arctan2(azCursor.pos[1]-azPlot.pos[1], 
                                azCursor.pos[0]-azPlot.pos[0])+np.pi/2
            elevation = (elPlot.bottom-elCursor.pos[1])/elPlot.height*140-50
            elevation = elevation*np.pi/180
            dist = (rCursor.pos[0]-rPlot.left)/rPlot.width*(maxR-minR) + minR
        else:
            if playing:
                scrubTime = time.time()-startTime-downTime
            if scrubTime >= limitTime:
                if posIn >= len(sourcePath):
                    break
                azimuth = sourcePath[posIn][0]
                elevation = sourcePath[posIn][1]
                dist = sourcePath[posIn][2]
                limitTime = limitTime + sourcePath[posIn][3]

                posIn += 1

        paz = -azimuth
        pel = elevation
        pr = dist

        # If moved bars, move radio
        if activeCursor in cursorList[1:3]:
            polCursor.pos = (int(np.cos(azimuth-np.pi/2)*polPlot.radius*(90-elevation*180/np.pi)/140)+polPlot.pos[0], 
                             int(np.sin(azimuth-np.pi/2)*polPlot.radius*(90-elevation*180/np.pi)/140)+polPlot.pos[1])

        ### DRAWS
        # if first:
        updateGraphics()
        # first=False
        clock.tick(FPS)


gameGUI(sourcePath)

pygame.quit()

stream.close()
wf.close()
p.terminate()

##### SAVE THE FILE
# fig = plt.figure()
# ax = fig.add_subplot(211)
# ax.plot(soundPlot[0])
# ax = fig.add_subplot(212)
# ax.plot(np.arange(len(soundPlot[1]))/44100, soundPlot[1])

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(311)
# ax2.plot(hrtfList[1][0])
# ax22 = fig2.add_subplot(312)
# ax22.plot(hrtfList[2][0])
# ax23 = fig2.add_subplot(313)
# ax23.plot(hrtfList[3][0])
# plt.show()
if recording:
    WAVE_OUTPUT_FILENAME = file[:-4] # filePath + fileName[:-4]
    WAVE_OUTPUT_FILENAME = WAVE_OUTPUT_FILENAME + ' hrtf.wav'
    # WAVE_OUTPUT_FILENAME = '/'.join(WAVE_OUTPUT_FILENAME)

    waveOut = wave.open(WAVE_OUTPUT_FILENAME, 'wb')
    waveOut.setnchannels(2)
    waveOut.setsampwidth(wf.getsampwidth())
    waveOut.setframerate(wf.getframerate())
    waveOut.writeframes(b''.join(recording))
    waveOut.close()
