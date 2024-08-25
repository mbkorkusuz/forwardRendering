# Rasterizer Render
## **Introduction**
This repository contains the implementation of well known forward rendering method in computer graphics.

## **Installation**
Make sure you have `imagemagick` on your system for image convertion.
If not visit [https://imagemagick.org/script/download.php]

1. Clone repository to your local machine
 ````text
git clone https://github.com/mbkorkusuz/forwardRendering.git
````
2. Navigate to the project directory
3. Compile the files
 ````text
make
```` 
5. Run `rasterizer` executable on input files, E.g
 ````text
./rasterizer inputs/sample.xml
````
    
## **Example Outputs**
Here are some rendered 3D objects
<div class="header">
  <h1>
    Box
  </h1>
</div>
<img src="/outputs/filled_box_4.ppm.png" alt="Box" title="Box" width=50% height=50%>
<div class="header">
  <h1>
    Horse and Mug
  </h1>
</div>
<img src="/outputs/horse_and_mug_2.ppm.png" alt="Horse and Mug" title="Horse and Mug" width=50% height=50%>
<div class="header">
  <h1>
    Flag
  </h1>
</div>
<img src="/outputs/flag_eu_1.ppm.png" alt = "Flag" title="Flag" width=50% height=50%>
