#different dimension in grf.c and other .c progs in img_estimator
path="/home/samir/simmulti/diffmulti/"
path1="/home/samir/simulation2/img_estimator/srcbin"
export path
export path1

#cp visnsedata_imgtge/vis1.03nse1.fits vis.FITS
cp vis_sc.fits vis.FITS

for i in `seq 1 1`;
do
    sed "s/-850/-2$i"76"/" input.grf >input.grfcp
    $path/grf input.grfcp tmp/img.FITS
  
    #cp visnsedata_imgtge/vis1.03nse$i.fits visnse.FITS
    #$path/visfitsgrid tmp/img.FITS visnse.FITS 1.
    #$path1/gridvisfits visnse.FITS tmp/tmp.FITS input.gridfits

    $path/visfitsgrid tmp/img.FITS vis.FITS 0.
    $path1/gridvisfits vis.FITS tmp/tmp.FITS input.gridfits

    $path1/collapse tmp/tmp.FITS tmp/tmpcollp.FITS
    $path1/image tmp/tmpcollp.FITS tmp/SI.FITS 7
    

    #theta_w is the radius in arcmin
    #$path1/window tmp/SI.FITS tmp/SW.fits 3 15.
    #$path1/rimage tmp/SW.fits tmp/BW.FITS 3
    #$path1/wt2 tmp/BW.FITS 4
    #$path1/image tmp/BW.FITS tmp/SW.fits 4

    $path1/multiply tmp/SI.FITS tmp/SW.fits tmp/SIW.FITS 7
    $path1/rimage tmp/SIW.FITS tmp/BIW.FITS 7
    
    #$path1/imgps tmp/BIW.FITS tmp/Mg$i.fits

    $path1/imgps tmp/BIW.FITS tmp/GVdiff$i.FITS
    $path1/binimgps tmp/GVdiff$i.FITS input.binimgps tmp/powerspec_diff$i.dat avgMg.fits

    #rm -rf tmp/*FITS
    #rm visnse.FITS
    
done

rm vis.FITS
#$path1/avgmg_img tmp/Mg avgMg.fits 128