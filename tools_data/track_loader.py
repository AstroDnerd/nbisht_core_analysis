from starter2 import *
reload(dl)
import looper2
reload(looper2)
import track_info
reload(track_info)


#Importing these files registers their contents in track_info.tracks
import track_files.tracks_u500
reload(track_files.tracks_u500)
import track_files.tracks_t200
reload(track_files.tracks_t200)
import track_files.tracks_t000
reload(track_files.tracks_t000)
import track_files.tracks_b000
reload(track_files.tracks_b000)
import track_files.tracks_u600
reload(track_files.tracks_u600)
import track_files.tracks_m000
reload(track_files.tracks_m000)
import track_files.tracks_u900
reload(track_files.tracks_u900)

if 'tracks' not in dir():
    tracks={}
    loops=tracks #for backwards compatibility

def load_tracks(namelist):
    if type(namelist) is not list:
        namelist = [namelist]
    for trackname in namelist:
        if trackname not in tracks:
            print("Load Loop",trackname)
            track=track_info.tracks[trackname]
            track_file = track.track_file
            sim_directory = track.sim_directory
            track = looper2.load_looper( track_file, directory=sim_directory, mode_fname=track.mode_fname)
            tracks[trackname]=track

