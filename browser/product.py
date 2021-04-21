
from starter2 import *
from collections import defaultdict

class plot():
    """container object that connects images with parameters"""
    def __init__(self,fname=".",parameters={}):
        self.fname=fname
        self.parameters=parameters

class product():
    """A tool to collect plots.  Takes
    myglob: a list of files. It's kind of redundant.
    regexp: a regular expression to match and extract parameters from.  For example, 
                /scratch/user/images/projection_x_c0123_n0100.png
            would want a regular expression to extract the core id (c0123) and frame (n0100)
    name: title on the column
    style: how to display the image.  Options are
        single: just make an img tag with each file name
        value: also then takes 
               fname: name of hdf5 file contating records
               field: name of field.  
               It expexts a record that looks like
                   fname.h5[field]['core_id']
            """
    def __init__(self, name="P", regexp="reg",myglob="glob",
                 parameters=['core_id','frame'],style='single',width=200,
                 fname=None,field=None):
        self.myglob=myglob
        self.regexp=regexp
        self.name=name
        self.parameters=parameters
        self.plots=defaultdict(lambda: list())
        self.style=style
        self.field=field
        self.fname=fname
        if style=='single':
            self.render = self.single_render
        elif style == 'core_id':
            self.render = self.core_id_render
        elif style == 'frames':
            self.render = self.frame_render
        elif style == 'value':
            self.render = self.value_render
        elif style == 'numbertest':
            self.render = self.number_render
        else:
            self.render = None
        
        self.width=width

    def render_head(self):
        return "<th> %s </th>"%self.name



    def check_glob(self):
        print("check glob")
        file_list = glob.glob(self.myglob)
        print(self.myglob)
        print(file_list)

    def get_frames(self):
        file_list = glob.glob(self.myglob)
        print(self.myglob)
        for fname in file_list:
            match = self.regexp.match(fname)
            if match is None:
                print("Error with regexp match")
                raise
            mygroups = match.groups()
            params = dict(zip(self.parameters,mygroups))
            core_id = int(params['core_id'])
            myplot = plot(fname,params)
            self.plots[core_id].append(myplot)
#       if len(self.parameters) > 1 and len(self.plots[core_id]) > 1:
#           for p in self.parameters:
#               if p != 'core_id':
#                   sort_key = p
#                   break
#           self.plots[core_id] = sorted( self.plots[core_id],key= lambda plot : plot.parameters[sort_key])

    def core_id_render(self,core_id):
        return "<td>%s</td>"%core_id
    def frame_render(self,core_id):
        if len( self.plots[core_id]) == 0:
            img = "x"
        else:
            if len(self.parameters) > 1 and len(self.plots[core_id]) > 1:
                for p in self.parameters:
                    #get a sort key that isn't core_id.
                    if p != 'core_id':
                        sort_key = p
                        break
                self.plots[core_id] = sorted( self.plots[core_id],key= lambda plot : plot.parameters[sort_key])
            img_tag_template = '<a h<figure><a href="%s"><img src="%s" width=%s id = %s></a><figcaption>%s</figcaption></figure>\n'
            fname = self.plots[core_id][0].fname
            caption = ""
            myid = "%s_c%04d"%(self.name, core_id)
            mynext = "%s_c%04d_next"%(self.name, core_id)
            myback = "%s_c%04d_back"%(self.name, core_id)
            img = img_tag_template%(fname, fname, self.width,myid,"")
            for nplot,plt in enumerate(self.plots[core_id]):
                #img += "<button onclick=set_image(%s,'%s')> n%04d</button>\n"%(myid,plt.fname,int(plt.parameters['frame']))
                next_frame = nplot+1
                if next_frame == len(self.plots[core_id]):
                    next_frame=0
                back_fname = self.plots[core_id][nplot-1].fname
                next_fname = self.plots[core_id][next_frame].fname
                img += "<button onclick=set_image('%s','%s')> n%04d</button>\n"%(myid,plt.fname,int(plt.parameters['frame']))
                print("render c%04d %s"%(core_id,str(plt.parameters)))
        out1 = "<td>%s</td>"%img
        return out1

    def value_render(self,core_id):
        if 'values' not in self.__dict__:
            fptr = h5py.File(self.fname,'r')
            values = fptr[self.field][()]
            core_ids = fptr['core_ids'][()]
            self.values = dict(zip(core_ids,values))
            fptr.close()
            #print(self.values)
        self.format='%0.2e'
        if core_id in self.values:
            out_str = self.format%self.values[core_id]
        else:
            out_str = 'xx'
        return "<td>%s </td>"%out_str


    def number_render(self,core_id):
        out = "<td>%0.2f</td>"%np.random.random()

        return out
    def single_render(self,core_id):
        #out1 = "<td>single render (%s)</td>"%self.style
        #img_tag = img_tag_template%(self.plots[core_id].fname)
        if len( self.plots[core_id]) == 0:
            img = "x"
        else:
            img_tag_template = '<a h<figure><a href="%s"><img src="%s" width=%s></a><figcaption>%s</figcaption></figure>'
            fname = self.plots[core_id][0].fname
            img = img_tag_template%(fname, fname, self.width,"")

        out1 = "<td>%s</td>"%img

        return out1 


