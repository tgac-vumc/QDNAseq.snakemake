#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cp -r $DIR/lb2 ./

echo "
<!DOCTYPE html>
<html lang="en-us">
<head>
  <meta charset="utf-8">
  <title>QDNAseq Profile Lightbox</title>
  <link rel="stylesheet" href="lb2/css/lightbox.css">
</head>
<body>

  <section>
    <div>
"

for i in *.png
do
echo "<a class=\"example-image-link\" href=\"${i}\" data-lightbox=\"profiles\"><img class=\"example-image\" src=\"${i}\" width="200" alt=\"\"/></a>"
done

echo "
</div>
  </section>

  <script src="lb2/js/lightbox-plus-jquery.min.js"></script>

</body>
</html>
"
