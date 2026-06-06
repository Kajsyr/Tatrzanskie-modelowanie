# Przewidywanie Występowania Wybranych Gatunków Alpejskich

## 1. Pozyskanie warstw z Google Earth Engine (GEE)

W celu modelowania niszy ekologicznej i przewidywania występowania wybranych gatunków alpejskich, zaimplementowano serię skryptów w środowisku Google Earth Engine (GEE). Pozwoliły one na wygenerowanie zestawu ciągłych zmiennych środowiskowych (predyktorów) w wysokiej rozdzielczości przestrzennej (10 m – 30 m), przyciętych do poligonu Alp.

Poniżej znajduje się szczegółowy opis pozyskanych komponentów środowiskowych wraz z kodem źródłowym:

### A. Wskaźniki Spektralne i Teledetekcyjne (Roślinność i Wilgotność)
Warstwy obliczone na podstawie danych satelitarnych **Sentinel-2 (Copernicus)** jako mediana z okresu letniego (01.06.2023 – 31.08.2023) dla pikseli o zachmurzeniu poniżej 20%. 

* **NDVI** (*Normalized Difference Vegetation Index*) – określa ogólną kondycję i gęstość biomasy.
* **EVI** (*Enhanced Vegetation Index*) & **SAVI** (*Soil-Adjusted Vegetation Index*) – wskaźniki zoptymalizowane pod kątem redukcji wpływu odsłoniętej gleby i efektów atmosferycznych, kluczowe w rzadkiej roślinności wysokogórskiej.
* **NDWI** (*Normalized Difference Water Index*) – odzwierciedla uwodnienie roślinności i wilgotność podłoża.
* **BSI** (*Bare Soil Index*) – identyfikuje obszary nagich skał, piargów i gleb pozbawionych szaty roślinnej.

<details>
<summary><b>Rozwiń kod GEE (Sentinel-2 Indices)</b></summary>

```javascript
//Shapefile alpy
var alpy_shp = ee.FeatureCollection("projects/lofty-complex-481222-k8/assets/alpy_shp");
var region = alpy_shp.geometry();

var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(region)
    .filterDate('2023-06-01', '2023-08-31')
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .median();

//Liczymy wskaźniki
//NDVI
var ndvi = s2.normalizedDifference(['B8', 'B4'])
  .clip(region)
  .rename('NDVI');
//EVI 
var evi = s2.expression('2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
  'NIR': s2.select('B8'), 'RED': s2.select('B4'), 'BLUE': s2.select('B2')
}).clip(region).rename('EVI');

//SAVI
var savi = s2.expression('((NIR - RED) / (NIR + RED + 0.5)) * 1.5', {
  'NIR': s2.select('B8'), 'RED': s2.select('B4')
}).clip(region).rename('SAVI');

//NDWI
var ndwi = s2.normalizedDifference(['B8', 'B11'])
  .clip(region)
  .rename('NDWI');

//BSI
var bsi = s2.expression('((SWIR + RED) - (NIR + BLUE)) / ((SWIR + RED) + (NIR + BLUE))', {
  'SWIR': s2.select('B11'), 'RED': s2.select('B4'),
  'NIR': s2.select('B8'), 'BLUE': s2.select('B2')
}).clip(region).rename('BSI');

// Centrowanie mapy do alp
Map.centerObject(region, 6);

// Definiowanie palet kolorów
var paleta_roslinnosc = ['#ffffff', '#ce7e45', '#dfc27d', '#aac98c', '#3e8f4c', '#1a491e'];
var paleta_wilgotnosc = ['#ece7f2', '#a6bddb', '#3690c0', '#023858'];
var paleta_skaly = ['#4d4d4d', '#bababa', '#e0e0e0', '#fddbc7', '#b2182b'];

// Dodawanie już przyciętych warstw do mapy
Map.addLayer(ndvi, {min: 0, max: 0.8, palette: paleta_roslinnosc}, 'NDVI');
Map.addLayer(evi, {min: 0, max: 1.0, palette: paleta_roslinnosc}, 'EVI');
Map.addLayer(savi, {min: 0, max: 0.8, palette: paleta_roslinnosc}, 'SAVI');
Map.addLayer(ndwi, {min: -0.3, max: 0.5, palette: paleta_wilgotnosc}, 'NDWI');
Map.addLayer(bsi, {min: -0.2, max: 0.4, palette: paleta_skaly}, 'BSI');

//Łączymy warstwy
var ecoFeatures = ee.Image.cat([
  ndvi, 
  evi, 
  savi, 
  ndwi, 
  bsi
]).float();

// Eksport całego pakietu zmiennych do jednego Assetu
Export.image.toAsset({
  image: ecoFeatures,
  description: 'eco_features_alpy_2023',
  assetId: 'eco_features_alpy_2023',
  scale: 10,
  region: region,
  maxPixels: 1e13
});
```
<details>

## Zaleganie śniegu
Czas zalegania okrywy śnieżnej ma fundamentalne znaczenie dla wegetacji roślin alpejskich. Warstwa przedstawia procentowy czas zalegania śniegu w okresie wiosenno-letnim (01.03.2023 – 30.06.2023). Skrypt maskuje chmury przy użyciu pasma SCL (Scene Classification Layer), oblicza wskaźnik NDSI (Normalized Difference Snow Index), a następnie wyznacza średnią częstość występowania pokrywy śnieżnej (NDSI > 0.4).

<details>
<summary><b>Rozwiń kod GEE (Sentinel-2 Indices)</b></summary>

```java
//Shapefile alpy
var alpy_shp = ee.FeatureCollection("projects/lofty-complex-481222-k8/assets/alpy_shp");
var region = alpy_shp.geometry();

// FUNKCJA MASKUJĄCA CHMURY I CIEŃ (Używa natywnego pasma SCL z Sentinel-2)
function maskCloudsS2(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  return image.updateMask(mask);
}

// Dane z copernicus usuwamy chmury
var s2_snow = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(region)
    .filterDate('2023-03-01', '2023-06-30')
    .map(maskCloudsS2); // Maskujemy pojedyncze chmury na każdym zdjęciu

// Funkcja obliczająca NDSI
var addNDSI = function(image) {
  var ndsi = image.normalizedDifference(['B3', 'B11']).rename('NDSI');
  return image.addBands(ndsi);
};

// Obliczenie prawdopodobieństwa wystąpienia śniegu (BEZ NaN)
var snowProbability = s2_snow.map(addNDSI)
    .map(function(img) {
      return img.select('NDSI').gt(0.4).rename('snow_mask');
    })
    .mean() 
    .multiply(100)
    .byte()
    .unmask(0) // ewentualne pozostałe braki danych (NaN) na wartość 0
    .rename('Snow_Duration_Pct');

//WIZUALIZACJA NA MAPIE
var snowVis = {
  min: 0,
  max: 100,
  palette: ['000000', '0d0887', '7e03a8', 'cc4778', 'f89441', 'eff821', 'ffffff']
};

Map.centerObject(alpy_shp, 6);
Map.addLayer(snowProbability.clip(region), snowVis, 'Czas zalegania śniegu (%) - Bez NaN');
```
## Właściwości Glebowe
Właściwości edaficzne pobrane z fizykochemicznej bazy danych SoilGrids (średnia dla warstwy powierzchniowej 0-5 cm) oraz globalnej mapy litologicznej USGS. Pozwalają one ująć w modelu preferencje podłoża u roślin (np. gatunki acydofilne i kalcyfilne).

* Soil_pH – odczyn gleby (wartości przemnożone przez 10).

* Soil_Organic_Carbon – zawartość węgla organicznego w glebie (indykator miąższości próchnicy).

* Soil_Sand - zawartość piasku
* Soil_Clay – zawartość gliny

<details>

<summary>Kod GEE<summary>

```java
var alpy_shp = ee.FeatureCollection("projects/lofty-complex-481222-k8/assets/alpy_shp");
var region = alpy_shp.geometry();

// pH (0-5cm) - wartości x10
var ph = ee.Image("projects/soilgrids-isric/phh2o_mean")
    .select('phh2o_0-5cm_mean')
    .clip(region)
    .rename('Soil_pH');

// Węgiel organiczny (0-5cm) - próchnica
var soc = ee.Image("projects/soilgrids-isric/soc_mean")
    .select('soc_0-5cm_mean')
    .clip(region)
    .rename('Soil_Organic_Carbon');

// Tekstura piasek
var sand = ee.Image("projects/soilgrids-isric/sand_mean")
    .select('sand_0-5cm_mean')
    .clip(region)
    .rename('Soil_Sand');

//Tekstrua glina
var clay = ee.Image("projects/soilgrids-isric/clay_mean")
    .select('clay_0-5cm_mean')
    .clip(region)
    .rename('Soil_Clay');

//Skały macierzyste
var usgs_litho = ee.Image("USGS/Global_Lithology/V1")
    .select('b1') 
    .clip(region)
    .rename('Lithology');

var soilLithoStack = ee.Image.cat([ph, soc, sand, clay, usgs_litho]).clip(region).float();

```
    
</details>


## Morfometria Terenu 

Na bazie numerycznego modelu wysokościowego SRTM (30 m) obliczono zestaw morfometrycznych pochodnych terenu, determinujących mikroklimat górski, grawitacyjne przemieszczanie mas oraz dystrybucję wilgoci:

* DEM – wysokość nad poziomem morza (kluczowy gradient termiczny).

* Slope – spadek terenu w stopniach.

* Aspect – ekspozycja stoku (kierunki świata w zakresie 0-360°).

* TRI (Terrain Ruggedness Index) – wskaźnik chropowatości terenu (identyfikacja m.in. ścian skalnych i stabilnych platform).

* Curvature – krzywizna terenu (różnicowanie form wklęsłych i wypukłych).

* TWI (Topographic Wetness Index) – wskaźnik potencjalnego uwodnienia topograficznego (akumulacja wód opadowych w dolinach i zagłębieniach).



