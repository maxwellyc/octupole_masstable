#==============================
#Defines the output to generate
#==============================
show loadpath
set terminal postscript eps enhanced color "bold" 35
set output "Broad_Grid_Data/Deformations/SLy4beta3.eps"   #Location/Name of the output file
#=====================================================
#Defines properties of the deformation values (z-axis)
#=====================================================
set cbrange [0:0.2]    #Range of the colorbox values
set cbtics 0,0.05,0.2  #Location of colorbox numerical labels
#======================================
#Defines properties of the x and y axes
#======================================
set xrange [0:304]     #Range of x values
set xtics 0,20,302     #Location of x-axis numerical labels
set yrange [0:124]     #Range of y values
set ytics 0,20,122     #Location of y-axis numbercal labels
set zrange [0:0.2]
#=============================================================
#Sets the size of the margins between the graph and the canvas
#=============================================================
set lmargin 0    #Size of the left margin between the graph and the canvas
set rmargin 0    #Size of the right margin between the graph and the canvas
set tmargin 0    #Size of the top margin between the graph and the canvas
set bmargin 0    #Size of the bottom margin between the graph and the canvas
#=======================
#Properties of the graph
#=======================
set size 2.975,1.8            #Size of the graph
set origin 1.00,0.15          #Sets the origin of the graph
set view map scale 1.0        #Displays the graph as a map
set xlabel "neutron number"   #Label of the x-axis
set ylabel "proton number"    #Label of the y-axis
#===========================================
#Properties of the colorbox (for the z-axis)
#===========================================
set palette defined (0 "white", 0.2 "red")
set colorbox user
set colorbox origin 3.0,0.55 size 0.1,0.5 
set cblabel "{/=60 {/Symbol b}_3}"
#================================
#Plot the beta_3 deformation data
#================================
#plot "Broad_Grid_Data/Ground_States/SLy4Axialz002to120.dat" using 2:1 notitle with points pointtype 5 pointsize 1.3,\
"Broad_Grid_Data/Ground_States/SLy4Axialz002to120.dat" using 2:1 notitle with points pointsize 1.5 pointtype 64 linecolor rgb "black" linewidth 3.0
splot "Broad_Grid_Data/Ground_States/SLy4Axialz002to120.dat" using 2:1:4 with points notitle pointtype 5 pointsize 1.3 palette linewidth 0.5,\
"Broad_Grid_Data/Ground_States/SLy4Axialz002to120.dat" using 2:1:4 notitle with points pointsize 1.5 pointtype 64 linecolor rgb "black" linewidth 1.0   #Outline of each color point
#####################################
#	UNEDF1
#
#set size 2.8,1.8
#set origin 1.00,1.45
#set view map
#
#
#set form x ""
#unset xlabel
#set ylabel "proton number
#unset colorbox
#
#
#set palette defined (-0.3 "blue", 0 "white", 0.4 "red")
#unset colorbox
#set label  "{/=60 UNEDF1}"  at 240,8
#set cblabel "{/=60 {/Symbol=b}_2}"
#

#sp "UNEDF1nuclei.dat" u 2:1:3  not  w p ps 1.3 pt 5 palette lw 0.5,\
#"UNEDF1nuclei.dat" u 2:1:3 not  w p ps 1.5 pt 64 lc rgb "black" lw 1.0

#unset label

#########################################
#	SLy4
#
#set size 2.8,1.8
#set origin 1.92,1.45
#set view map
#
#
#set form x ""
#set form y ""
#unset xlabel
#unset ylabel
#unset colorbox
#
#set palette defined (-0.4 "blue", 0 "white", 0.4 "red")
#unset colorbox
#
#unset colorbox
#set label  "{/=60 SLy4}"  at 240,8


#sp "Sep.SLY4-LN.dat" u 4:5:($5<122&$127>0&$128>0?$59-$58:1/0)  not  w p ps 1.3 pt 5 palette lw 0.5,\
#"Sep.SLY4-LN.dat" u 4:5:($5<122&$127>0&$128>0?1:1/0) not  w p ps 1.5 pt 64 lc rgb "black" lw 1.0


#unset label

###########################################
#	SV-min
#
#set size 2.8,1.8
#set origin 1.92,0.15
#set view map
#
#
#set form x 
#unset ylabel
#set xlabel "neutron number"
#
#
#set palette defined (-0.3 "blue", 0 "white", 0.4 "red")
#unset colorbox
#set label  "{/=60 SV-min}"  at 240,8
#set cblabel "{/Symbol=b}_2"
#
#sp "Sep.SV-MIN-LN.dat" u 4:5:($5<122&$127>0&$128>0&$59<0?$59-$58:1/0)  not  w p ps 1.3 pt 5 palette lw 0.5,\
#"Sep.SV-MIN-LN.dat" u 4:5:($5<122&$127>0&$128>0&$59<0?1:1/0) not  w p ps 1.5 pt 64 lc rgb "black" lw 1.0


#######################

#unset label

###########################
#	UNEDF2
#
#set size 2.8,1.8
#set origin 1.00,2.75
#set view map
#
#
#set form x ""
#set form y
#unset xlabel
#unset ylabel
#set ylabel "proton number
#unset colorbox
#
#
#set palette defined (-0.3 "blue", 0 "white", 0.4 "red")
#unset colorbox
#set label  "{/=60{UNEDF2}}"  at 240,8
#set cbrange [-0.2:0.6]
#

#sp "UNEDF2nuclei.dat" u 2:1:3  not  w p ps 1.3 pt 5 palette lw 0.5,\
#"UNEDF2nuclei.dat" u 2:1:3 not  w p ps 1.5 pt 64 lc rgb "black" lw 1.0


#unset label

###########################
#	SkP
#
#set size 2.8,1.8
#set origin 1.92,2.75
#set view map
#
#
#set form x ""
#set form y ""
#unset xlabel
#unset ylabel
#
#set palette defined (-0.4 "blue", 0 "white", 0.4 "red")

#
#sp "Sep.SKP-LN.dat" u 4:5:($5<122&$127>0&$128>0?$59-$58:1/0)  not  w p ps 1.3 pt 5 palette lw 0.5,\
#"Sep.SKP-LN.dat" u 4:5:($5<122&$127>0&$128>0?1:1/0) not  w p ps 1.5 pt 64 lc rgb "black" lw 1.0




#####################
