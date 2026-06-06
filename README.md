# Przewidywanie Występowania Wybranych Gatunków Alpejskich

## Pozyskanie warstw z Google Earth Engine (GEE)

W celu modelowania niszy ekologicznej i przewidywania występowania wybranych gatunków alpejskich, zaimplementowano serię skryptów w środowisku Google Earth Engine (GEE).
Poniżej znajduje się szczegółowy opis pozyskanych komponentów środowiskowych wraz z kodem źródłowym:

### Wskaźniki Spektralne i Teledetekcyjne (Roślinność i Wilgotność)
Warstwy obliczone na podstawie danych satelitarnych **Sentinel-2 (Copernicus)** jako mediana z okresu letniego (01.06.2023 – 31.08.2023) dla pikseli o zachmurzeniu poniżej 20%. 

* **NDVI** (*Normalized Difference Vegetation Index*) – klasyczny wskaźnik ilości biomasy.
  <center>$$NDVI = \frac{NIR - RED}{NIR + RED}$$</center>
  
* **EVI** (*Enhanced Vegetation Index*) – wskaźnik zoptymalizowany pod kątem redukcji wpływu efektów atmosferycznych i nasycenia wysoką biomasą.
  <center>$$EVI = 2.5 \times \frac{NIR - Red}{NIR + 6 \times Red - 7.5 \times Blue + 1}$$</center>

* **SAVI** (*Soil-Adjusted Vegetation Index*) – wskaźnik minimalizujący wpływ jasności tła gleby w obszarach o rzadkiej roślinności.
  <center>$$SAVI = \frac{NIR - Red}{NIR + Red + L} \times (1 + L)$$</center>
  
* **NDWI** (*Normalized Difference Water Index*) – odzwierciedla uwodnienie roślinności i wilgotność podłoża.
  <center>$$NDWI = \frac{Green - NIR}{Green + NIR}$$</center>
  
* **BSI** (*Bare Soil Index*) – identyfikuje obszary nagich skał, piargów i gleb pozbawionych szaty roślinnej.
  <center>$$BSI = \frac{(SWIR + Red) - (NIR + Blue)}{(SWIR + Red) + (NIR + Blue)}$$</center>


<details>
<summary> Kod GEE </summary>

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
</details>

### Zaleganie śniegu
Czas zalegania śniegu. Okres wiosna-lato z 2023. Wynik w procentach. 0 to brak śniegu. Bazujemy na wskaźniku NDSI (Normalized Difference Snow Index). 
<center>$$NDSI=\frac{GREEN-SWIR}{GREEN+SWIR}$$</center>

<details>
<summary>Kod GEE</summary>

```javascript
//Shapefile alpy
var alpy_shp = ee.FeatureCollection("projects/lofty-complex-481222-k8/assets/alpy_shp");
var region = alpy_shp.geometry();

// FUNKCJA MASKUJĄCA CHMURY I CIEŃ 
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

// Obliczenie prawdopodobieństwa wystąpienia śniegu
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
Map.addLayer(snowProbability.clip(region), snowVis, 'Czas zalegania śniegu (%));
```
</details>

### Właściwości Glebowe
Właściwości edaficzne pobrane z fizykochemicznej bazy danych SoilGrids (średnia dla warstwy powierzchniowej 0-5 cm) oraz globalnej mapy litologicznej USGS. Pozwalają one ująć w modelu preferencje podłoża u roślin (np. gatunki acydofilne i kalcyfilne).

* Soil_pH – odczyn gleby (wartości przemnożone przez 10).

* Soil_Organic_Carbon – zawartość węgla organicznego w glebie (indykator miąższości próchnicy).

* Soil_Sand - zawartość piasku
  
* Soil_Clay – zawartość gliny

<details>
<summary>Kod GEE</summary>

```javascript
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

### Skład Mineralny
Ze względu na braki w danych geomorfologicznynych używamy równeiż danych z satelity które zapewniąja pełniejszy obraz. Szczególnie zależey nam na rozróznieniu pomiędzy granitami a wapieniami.

* **Silicates** – wskaźnik krzemianów, charakterystyczny dla granitów.

* **Carbonates** – wskaźnik węglanów, charakterystyczny dla wapieni

* **Iron_Oxides** – wskaźnik tlenków żelaza

<details>
<summary> GEE KOD </summary>
    
```javascript
var alpy_shp = ee.FeatureCollection("projects/lofty-complex-481222-k8/assets/alpy_shp");
var region = alpy_shp.geometry();

var landsat_collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(region)
  .filterDate('2021-01-01', '2025-12-31') 
  .filter(ee.Filter.calendarRange(6, 9, 'month')) 
  .filter(ee.Filter.lt('CLOUD_COVER', 15));

var image = landsat_collection.median().clip(region);

var silicates = image.normalizedDifference(['SR_B6', 'SR_B7']).rename('Silicates').clip(region); 
var carbonates = image.normalizedDifference(['SR_B5', 'SR_B6']).rename('Carbonates').clip(region);
var ironOxides = image.normalizedDifference(['SR_B4', 'SR_B3']).rename('Iron_Oxides').clip(region); 

var highGeoStack = ee.Image.cat([silicates, carbonates, ironOxides]).float();
```
</details>

### Nasłonecznienie
Dane z satelity z poprawką na kształt terenu w oparciu o ekpozycje i skośność.

<details>
<summary>Kod GEE</summary>

```javascript
var alpy_shp = ee.FeatureCollection("projects/lofty-complex-481222-k8/assets/alpy_shp");
var region = alpy_shp.geometry();

var solar_collection = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_BY_HOUR")
    .filterBounds(region)
    .filterDate('2023-05-01', '2023-09-30')
    .select('surface_solar_radiation_downwards_hourly');

var totalSolar = solar_collection.mean().divide(1000000).clip(region);

var dem = ee.Image("USGS/SRTMGL1_003").clip(region);
var terrain = ee.Terrain.products(dem);
var slopeRad = terrain.select('slope').multiply(Math.PI / 180);
var aspectRad = terrain.select('aspect').multiply(Math.PI / 180);

var hli = dem.expression(
  'cos(aspect - 3.927) * sin(slope)', {
    'aspect': aspectRad, 'slope': slopeRad
  }
);

var solarModifier = hli.add(1).multiply(0.5).add(0.5);
var topo_solar = totalSolar.resample('bicubic').reproject({crs: dem.projection(), scale: 30}).multiply(solarModifier).rename('Topo_Solar_Radiation_MJ');
```
</details>


### Morfometria Terenu 

Na bazie numerycznego modelu wysokościowego SRTM (30 m) obliczono zestaw morfometrycznych pochodnych terenu, determinujących mikroklimat górski, grawitacyjne przemieszczanie mas oraz dystrybucję wilgoci:

* **DEM** – wysokość nad poziomem morza

* **Slope** – spadek terenu w stopniach.

* **Aspect** – ekspozycja stoku (kierunki świata w zakresie 0-360°).

* **TRI** (Terrain Ruggedness Index) – wskaźnik chropowatości terenu

* **Curvature** – krzywizna terenu różnicowanie form wklęsłych i wypukłych.

* **TWI** (Topographic Wetness Index) – wskaźnik potencjalnego uwodnienia topograficznego czyli gdzie spływa woda

<details>
<summary>Kod GEE</summary>
    
```javascript
var alpy_shp = ee.FeatureCollection("projects/lofty-complex-481222-k8/assets/alpy_shp");
var region = alpy_shp.geometry();

var dem = ee.Image("USGS/SRTMGL1_003").clip(region);
var terrain = ee.Terrain.products(dem);

var slope = terrain.select('slope').rename('Slope');
var aspect = terrain.select('aspect').rename('Aspect');
var tri = dem.convolve(ee.Kernel.fixed(3, 3, [
  [-1/8, -1/8, -1/8],
  [-1/8,  1,   -1/8],
  [-1/8, -1/8, -1/8]
])).abs().rename('TRI');

var curvature = ee.Terrain.slope(ee.Terrain.slope(dem)).rename('Curvature');
var slopeRad = slope.multiply(Math.PI).divide(180);
var twi = dem.focal_mean(300, 'circle', 'meters').divide(slopeRad.add(0.01).tan()).log().rename('TWI');

var topographyStack = ee.Image.cat([dem.rename('DEM'), slope, aspect, tri, curvature, twi]).float();
```
</details>



## Pozyskanie stanowisk roślin z bazy GBIF ()
Proces pobierania i przygotowania danych o lokalizacji gatunków zozstał zautomatyzowany przy pomocy potoku Nextflow. Składa się z następujących kroków:

Ekstrakcja nazw taksonów: Skrypt pobiera unikalną listę gatunków z arkusza Excel (.xlsx), a pipeline Nextflow rozdziela ją na niezależne, równoległe wątki.

### Pobranie z GBIF
Wybrane gatunki zostały pobrane używajac API GBIF (rgbif) w nastepujacych kryteriach:

* Rok obserwacji: od 1990 roku wzwyż.
* Niepewność współrzędnych: poniżej 20 metrów.
* Brak znanych błędów geoprzestrzennych i wymagana obecność współrzędnych GPS.
* Zgrubne ograniczenie do obszaru Europy Środkowej (ramka przestrzenna WKT).

### Filtrowanie
Nastepnie wyniki zostały poddane walidacji i zostały usunięte punkyt:
* z wartościami zerowymi (0,0)
* punkty pokrywające się z centrami miast i stolic,
* punkty w obszarach zurbanizowanych
* w lokalizacji ogrodów botanicznych, instytucji naukowych i herbariów
* poza obszarem Alp

Pobrano 17215 rekordów dla 26 gatunków. Po odfiltrowaniu uzyskano 12524 rekordów któ©e zostały przekazaen do dalszej analizy.

## Feature Selection
Następnie dla stanowisk określono listę cech z warstw z GEE oraz z Bioclim. Bazując na tych danych wybrano cechy które zostały użyte w dalszej części do tworzenie modeli. Szczegóły w [PLIKU](01_feature_selection.ipynb)
   
![Heatmap](plots/correlation_heatmap_features.png)
![Histogramy](plots/feature_histograms.png)
### Feature selesction East-West
W trakcie wyboru warstw przy pewnym nie radykalnym wyborze zawierajacym sporo warst klimatycznych udało nam się zauwazyć dość wyraźny podziął na oba klastey. Po sklastrowaniu ozkaząło sie że skorelowane są one z rozkładem geogreaficznym. Co ciekawe pokrywa się to z umowną granicą geograficzną pomiedy Alpami Wschodnimi i Zachodnimi wynikającą z geomorfologii.
![PCA](plots/pca_east_west.png)
![MAP](plots/mapa_east_west.png)
Po odcięciu warst skorelowanych ze sobą ten podział znika. Pytanie czy ten podział wynika ze wspóliniowosci i jest fałszywy czy porzzez zbyt rygorytyczne odcinanie warst tracimy informację biologiczną? Nie wiemy.

## Uczenie
Szczegóły etapu uczenia w pliku [02_Learning](02_Learning.ipynb)
Do uczenia wybrano tylko te gatunki, dla których mieliśmy więcej niż 500 obserwacji. Końcowa liczba gatunków używana do uczenia wyniosła 9 i zbiór ten był całkiem niezbalansowany (poniżej wykres):
![n_obserwacji](plots/hist_obserwacji.png).

### Przygotowanie rekordów do uczenia
Jako rekordów nie traktowano pojedynczych obserwacji.
Zamiast tego staraliśmy się zdefiniować "obszary". Dyskutowaliśmy dwie opcje: 1. Grid na mapie. Wszystkie obserwacje w gridzie zostają włączone jako obserwacje z grida. 2. DBscan po koordynatach.
Ostatecznie wybraliśmy opcję dbscan jako sensowniejszą geograficznie, choć ryzykowniejszą, gdyby obserwacje układały się w jakieś duże struktury. Natomiast liczba obserwacji w naszych klastrach była koniec końców raczej sensowna (<16). Liczba klastrów wyniosła 6905.
Zmienne środowiskowe dla klastra ustanowiliśmy jako średnią z obserwacji w klastrze. Dla większości naszych zmiennych jest to raczej sensowne przybliżenie klimatu w małych klastrach. Natomiast dla niektórych zmiennych, zwłaszcza tych związanych z ukształtowaniem terenu może gubić sensowne informacje.

W ten sposób dla wszystkich naszych rekordów istnieje przynajmniej jeden gatunek, który został "w tym rekordzie" zaobserwowany.
Mamy, więc pewność, że był tam jakiś obserwator. Jeżeli jednocześnie nie mamy z tego miejsca rekordu o wystąpieniu innego gatunku, możemy przyjąć, że go tam nie ma. Oczywiście w ten sposób traktujemy to jako układ statyczny, bez żadnej strzałki czasu.
Drugie założenie, które w ten sposób poczyniamy to przyjęcie braku interakcji między gatunkami.
Można to założenie zaatakować sprawdzając czy gatunki zajmują podobne nisze klimatyczne (czy dobrze się klastrują w tej domenie), a czy mimo to nie występują razem. Wtedy mamy np. znaczącą wskazówkę co do konkurencji między nimi.
Z drugiej strony moglibyśmy też po prostu uwzględnić obecność innych gatunków w danym miejscu jako zmienną do predykcji dla tego jednego gatunku.
Niestety nie zajmowaliśmy się tym tutaj, a jako robocze przybliżenie założyliśmy, że tych interakcji nie ma, choć spojrzeliśmy sobie na kilka metryk (korelację Pearsona, pointwise-mutual-information i prawdopodobieństwo warunkowe) z czystej ciekawości (heatmapy w 02_Learning.ipynb), ale bez uzupełnienia tego o klastrowanie po klimacie nie wiedzieliśmy co możemy z tymi wielkościami zrobić, a na porządną analizę tego typu zależności brakło nam czasu. Na pewno byłoby to ciekawe rozszerzenie projektu.

### Modele
Model trenowany na 5-foldach. Nie wydzielaliśmy danych testowych, bo uznaliśmy, że cross-validation jest wystarczająca. Choć możliwe, że jeden testowy (taki 100% agnostyczny) przydałby się do porównania między modelami.
Do trenowania używaliśmy RandomizedSearchCV z sklearn.

Wytrenowaliśmy trzy modele:
1. Random Forest, który przewiduje wektor z Y^N. Gdzie N to liczba gatunków.
2. Model RF podobny do 1, ale taki, w którym obszary z foldów testowych są z innych bloków (większych obszarów) niż obszary z foldów treningowych. Chcieliśmy w ten sposób uniknąć sytuacji, w której klaster 'z jednego końca polany' jest w secie treningowym, a klaster 'z drugiego końca tej samej polany' jest w secie testowym.
3. "N*Y" - Model, który tak naprawdę jest zlepkiem N modeli RF, w którym każdy model przewiduje jedną liczbę. Chcieliśmy porównać, jak to wpłynie na czas treningu oraz na jakość predykcji. Znowu groupkfold po blokach.



### Metryki
W czasie pierwszych próbnych podejść dostawaliśmy błąd na uczeniu spowodowany obliczaniem recall, albo precyzji. Doszliśmy do wniosku (nie bez pytania llma, że może chodzić o zera w testowym foldzie). Metryka w tego typu zadaniach jest nieoczywista. f1-score się nam nie spodobał. O tyle o ile precyzja TP/(TP+FP) jest całkiem sensowna, bo rośnie gdy mamy mało false positiwów, to recall nie jest aż tak interesujący bo mówi nam o propocji TP/(TP+FN). Duża liczba false negatiwów nas aż tak nie boli, jak mogłaby nas boleć duża liczba false positiwów.

Jako metrykę, którą optymalizujemy w modelu (argument scoring) wybraliśmy średnią dla wszystkich gatunków z wartości average_precision_score (z sklearn) uśredniony po wszystkich foldach, ale tych, w których zaobserwowaliśmy przynajmniej jedną jedynkę dla gatunku. Dla szczegółów odsyłamy do funkcji map_scorer().

average_precision_score (AP) to suma po wartościach precyzji dla różnych progów recallu. Suma ważona przez przyrost w recall. Czyli im więcej jedynek będzie miało wyższe precision niż zera tym AP będzie większe. Dla niezbalansowanych danych losowe rozrzucenie danych odpowiada AP=propocji jedynek. Dlatego będziemy normalizować:
lift = (AP-base)/(1-base), gdzie base to propocja jedynek. Dla losowego rozrzucenia mamy lift=0. Dla idealnego AP mamy lift=1.


### Feature importance dla pierwszego modelu
Ważność cech sprawdziliśmy w oparciu o metodę permutacyjną. Patrzyliśmy jak zmienia się wartość map_scorer() kiedy spermutujemy wartości danej cechy.
Ważność sprawdziliśmy jedynie dla pierwszego modelu (poniżej barplot).

![feat_importance](plots/feature_importance_rf1.png)

Jako najważniejsza cecha wybiła się bio_9, czyli "Mean Temp Driest Quarter". Co jest ciekawe samo w sobie, natomiast wynika to też z odrzucenia innych zmiennych podobnego typu w czasie feature selection.

Totalnie nieważne okazały się zmienne Slope i Curvature. Co ciekawe "wizualna" ocena na density plot wskazuje na podobieństwo rozkładów losowych punktów i punktów z obserwacji. Więc, być może to już samo w sobie było przesłanką do odrzucenia ich z uczenia. Z drugiej strony Slope i Curvature używane do uczenia to tak naprawdę średnia dla klastra, więc może, dla klastrów gdzie n>1, zatarła się ważna informacja, w co raczej wątpimy.

### Wyniki i wnioski

---

Suche wyniki wyprintowane w outputach notebooku 02_Learning.ipynb.
![dumbbell](plots//dumbbell.png)

Pierwsze co widać, to to, że pierwszy model dla większości gatunków wystrzelił wyższy wynik.
Drugi i trzeci model mają w miarę podobne wyniki, dlatego będziemy mówić tylko o pierwszym i drugim. Ale najpierw jeszcze zauważymy, że jest to dobra informacja, ponieważ oba modele Y^N trenowały się dużo krócej niż N*Y. Jest to ciekawy wynik do zapamiętania, jeżeli chodzi o optymalizację kodu wynik/czas. Choć pewnie nie będzie to prawda dla każdego zestawu danych. Tutaj się sprawdza.

Druga sprawa, która jest w naszej ocenie najciekawszym wynikiem to to, że lift spada znacząco kiedy pogrupujemy foldy tak, aby nie testować na klastrach, które są blisko klastrów treningowych. Czyli random-split może nam przekłamać model. Mało przekłamuje dla Linaria alpina co początkowo uznaliśmy za skutek tego, że jest to najliczniejszy gatunek, ale raczej nie jest to trend, więc nie będziemy takich wniosków tutaj wysuwać.
Chętnie pociągnęlibyśmy ten wątek dalej i patrzyli jak zmienia się ten spadek w zależności od przyjętych rozmiarów bloków (a być może nawet od przyjętego rozmiaru klastrowania).

---

Chcąc sprawdzić jakoś użyteczność tych modeli wrzuciliśmy do pierwszego z nich losowe punkty z mapy (wcześniej odfiltrowaliśmy te, które są blisko obserwacji używając BallTree).
Poniżej wykres liczby 1 według modelu dla różnego progu prawdopobieństwa:

![threshold](plots/threshold_on_random.png)

Widzimy coś bardzo informatywnego: Linaria alpina ma najmniej selektywne predykcje. Być może przez to, że ma wysoką obecność w danych wymaga dużego progu. Pozostałe modele są bardzo ostrożne i niskie progi już odcinają nam dużo "zer".
Podejrzewamy, że możnaby to jakoś ustabilizować dla Linaria używając tła z dużą liczbą zer lub obserwacji z innych gatunków, równie licznych lub takich, których byśmy nie modelowali, a jedynie brali punkty z mapy jako miejsca, w których był obserwator.
---


### Dywagacje
Być może różnice w spadkach między gatunkami da się wytłumaczyć jakimiś zmiennymi klimatycznymi. Ale nie sprawdziliśmy tego.

Ciężko nam oceniać nasze wyniki jako wyniki z modelowania niszy ekologicznej, bo nie znamy się na tym i brakowało nam porównania ze state-of-art modelami.

# Uwagi:
- Deduplikacja rekordów powinna nastąpić przed feature selection.
- Obserwacje zostały zbinaryzowane, deduplikacja + łączenie sąsiednich rekordów z tego samego gatunku. Zamiast tego
możnaby przypisać temu jakąś wagę. Jeżeli mamy więcej rekordów z tego miejsca to być może panują tam lepsze warunki. Z drugiej strony unikamy problemu powielonych obserwacji lub celowo pominiętych osobników, kiedy obserwator tak naprawdę spotkał większą ich liczbę.




