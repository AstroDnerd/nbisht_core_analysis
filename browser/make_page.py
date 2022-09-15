from starter2 import *
import product
reload(product)

head_dumb="""
<html>
<head>
<script src="./sorttable.js"></script>
<script>
function set_image(image_id,fname){
    this_image =     document.getElementById(image_id)
    this_image.src = fname;
}
</script>
 <link rel="stylesheet" href="sorttable.css">
</head>
<body>
"""

def make_page(product_list,core_list=None, htmlname="./output.html", captions=None):
    if core_list is None:
        core_list=[]
        for p in product_list:
            mycores=list(p.plots.keys())
            core_list+=mycores
        core_list=np.unique(np.array(sorted(core_list)))

    fptr = open(htmlname,'w')
    fptr.write(head_dumb)
    #loader=jinja2.FileSystemLoader('.')
    #env = jinja2.Environment(loader=loader)
    #main_template  = env.get_template('budget_template_1.html')
    fptr.write("Stuff<br>\n")
    fptr.write('<table class="sortable" id="example">\n')
    fptr.write('<tr>')
    fptr.write('<td>Core</td>')
    for p in product_list:
        fptr.write(p.render_head())
    fptr.write('</tr>\n')

    for nc,core_id in enumerate(core_list):
        fptr.write('<tr>')
        fptr.write('<th> %d </th>'%core_id)
        for p in product_list:
            fptr.write(p.render(core_id))
        fptr.write('</tr>\n')
    
    if captions is not None:
        fptr.write('<tr>')
        fptr.write('<th> Captions </th>')
        for nprod,p in enumerate(product_list):
            fptr.write("<td>")
            fptr.write( captions[nprod])
            fptr.write("</td>")
        fptr.write('</tr>\n')
                
    fptr.write("</table>\n")
    fptr.write("</body>\n")
    fptr.close()
    


def make_page_multi(product_list,core_list=None, htmlname="./output.html"):
    if core_list is None:
        core_list=[]
        for p in product_list:
            mycores=list(p.plots.keys())
            core_list+=mycores
        core_list=np.unique(np.array(sorted(core_list)))

    fptr = open(htmlname,'w')
    fptr.write(head_dumb)
    #loader=jinja2.FileSystemLoader('.')
    #env = jinja2.Environment(loader=loader)
    #main_template  = env.get_template('budget_template_1.html')
    fptr.write("Stuff<br>\n")
    fptr.write('<table class="sortable" id="example">\n')
    fptr.write('<tr>')
    fptr.write('<td>Core</td>')
    for p in list(product_list.values())[0]:
        fptr.write(p.render_head())
    fptr.write('</tr>\n')

    for simname in sorted( product_list.keys()):
        for nc,core_id in enumerate(core_list[simname]):
            fptr.write('<tr>')
            fptr.write('<th> %d </th>'%core_id)
            for p in product_list[simname]:
                fptr.write(p.render(core_id))
            fptr.write('</tr>\n')
    
                
    fptr.write("</table>\n")
    fptr.write("</body>\n")
    fptr.close()
    


