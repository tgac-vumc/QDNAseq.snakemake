#!/bin/bash

echo "
<!DOCTYPE html>
<html lang="en-us">
<head>
  <meta charset="utf-8">
  <title>QDNAseq Profile Lightbox</title>
  <link rel="stylesheet" href="http://ccagc-gen01.vumc.nl/js/lightbox2-master/dist/css/lightbox.css">
</head>
<body>

  <section>
    <div>
"

for i in `find ./ -iname "*.png"`
do
echo "<a class=\"example-image-link\" href=\"${i}\" data-lightbox=\"profiles\"><img class=\"example-image\" src=\"${i}\" width="200" alt=\"\"/></a>"
done

echo "
</div>
  </section>

  <script src="http://ccagc-gen01.vumc.nl/js/lightbox2-master/dist/js/lightbox-plus-jquery.min.js"></script>

</body>
</html>
"
