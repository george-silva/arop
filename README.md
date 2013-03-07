arop
====

Automated Registration and Orthorectification Package

*Note: this is a software that was developed at USDA (United States Department of Agriculture)*.
I'm just hosting the source to make it easier for others to consume it directly.
I'm not the author of this software.

Here is the link for the original [entry](http://globalchange.nasa.gov/KeywordSearch/Metadata.do?Portal=GCMD_Services&KeywordPath=ServiceParameters%7CDATA+ANALYSIS+AND+VISUALIZATION&EntryId=USDA_ARS_AROP&MetadataView=Full&MetadataType=1&lbnode=mdlb5)

---

## Abstract:

Accurate geo-referencing information is a basic requirement for combining remote satellite imagery with other geographic information. To detect changes in time-series satellite images, it is extremely important for the images to be precisely co-registered and orthorectified, so that images acquired from different sensors and dates can be compared directly. 

Precise registration relates satellite images to the ground reference based on carefully selected ground control points between the image and corresponding ground objects. Co-registration matches two images based on the tie points in the images. Although most commercial remote sensing software packages provide precise registration or co-registration, the selection of control points or tie points needs to be done manually. This is a time-consuming process that is not feasible for a large amount of data processing. Even though some satellite data are precisely registered and orthorectified, mismatches of up to a few pixels are common in early data sources. 

The topographical variations of the earth’s surface and the satellite view zenith angle affect the pixel’s distance projected onto the satellite image. The distortion inherent in the image is determined by topographical elevation. For example, an off-nadir view of 5 degrees at a site with a terrain height of 1000 meters can cause displacement of 120 meters (about 4 Landsat Thematic Mapper (TM) pixels). The orthorectification process is used to correct the pixel displacement caused by the topographical variations at the off-nadir viewing and to make the image orthographic, with every pixel in its correct location regardless of elevation and viewing direction. 

We developed an automated registration and orthorectification package (AROP) for Landsat and Landsat-like data processing. It uses precisely registered and orthorectified Landsat data (e.g., GeoCover or recently released free Landsat Level 1T data from the USGS EROS data center) as the base image to co-register and orthorectify the warp image and thus to make geo-referenced time-series images consistent in the geographic extent, spatial resolution, and projection. This package has been tested on the Landsat Multi-spectral Scanner (MSS), TM and Enhanced TM Plus (ETM+), Terra ASTER, CBERS CCD, and IRS-P6 AWiFS data. 

---

The development of the AROP package was supported by the U.S. Geological Survey (USGS) Landsat Science Team project and the NASA EOS project. The package was initially developed at the NASA Goddard Space Flight Center (from September 2005 to June 2011). Previous version can be downloaded from the NASA GSFC Landsat Ecosystem Disturbance Adaptive Processing System (LEDAPS) web site. Further improvement and continuous maintenance are now being undertaken in the Hydrology and Remote Sensing Laboratory, Agricultural Research Service, U.S. Department of Agriculture (USDA) by Dr. Feng Gao. 

