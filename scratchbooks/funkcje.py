import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler



# =====================================================================
# FUNKCJA 1: OBLICZENIA I DOPASOWANIE PCA
# =====================================================================
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def calculate_pca(
    df, target_col="true_present", n_components=9
):
    """Przeprowadza PCA na cechach numerycznych i łączy wyniki

    z oryginalnymi kolumnami nienumerycznymi.
    """
    # 1. Usuwamy wiersze z brakami danych
    df_clean = df.dropna()

    # 2. Separujemy klasę docelową oraz kolumny nienumeryczne
    y_clean = df_clean[target_col]

    # Pobieramy wszystkie kolumny nienumeryczne (oprócz targetu, jeśli nim był)
    df_non_numeric = df_clean.select_dtypes(exclude=[np.number]).drop(
        columns=[target_col], errors="ignore"
    )

    # Wybieramy tylko numeryczne cechy do PCA
    X_numeric = df_clean.select_dtypes(include=[np.number]).drop(
        columns=[target_col], errors="ignore"
    )

    # 3. Standaryzacja cech
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_numeric)

    # 4. Dynamiczne dopasowanie liczby komponentów
    actual_components = min(
        n_components, X_numeric.shape[0], X_numeric.shape[1]
    )

    # 5. Dopasowanie PCA
    pca = PCA(n_components=actual_components)
    X_pca = pca.fit_transform(X_scaled)

    # 6. Procent wyjaśnianej wariancji
    explained_variance = pca.explained_variance_ratio_ * 100

    # 7. Budowanie DataFrame z komponentami PC
    pc_columns = [f"PC{i+1}" for i in range(actual_components)]
    df_pca_results = pd.DataFrame(
        X_pca, columns=pc_columns, index=df_clean.index
    )

    # 8. Łączenie: nienumeryczne + komponenty PCA
    df_final = pd.concat([df_non_numeric, df_pca_results], axis=1)

    # Mapowanie klasy docelowej
    df_final["Target_Class"] = y_clean.map({0: "Pseudo-absence", 1: "Presence"})

    return df_final, explained_variance


# =====================================================================
# FUNKCJA 2: WIZUALIZACJA PRZESTRZENI PC1 I PC2
# =====================================================================
def plot_pca_space(
    df_pca,
    explained_variance,
    hue_col="Target_Class",  # Dynamiczna kolumna do kolorowania punktów
    title="Wizualizacja przestrzeni cech przy użyciu PC1 i PC2",
    palette=None,
    alpha=0.2,
    s=40,
):
    """Generuje wykres rozrzutu dla dwóch pierwszych składowych głównych (PC1 i PC2)."""

    # Sprawdzenie czy w przekazanym DataFrame istnieją wymagane składowe
    if "PC1" not in df_pca.columns or "PC2" not in df_pca.columns:
        print(
            "Błąd: Przekazany DataFrame nie zawiera składowych PC1 lub PC2. Uruchom najpierw calculate_pca."
        )
        return

    # Sprawdzenie czy wskazana kolumna do kolorowania istnieje w DataFrame
    if hue_col not in df_pca.columns:
        print(f"Błąd: Kolumna '{hue_col}' nie istnieje w przekazanym DataFrame.")
        return

    # Jeśli użytkownik nie podał własnej palety i koloruje domyślną klasę docelową
    if palette is None and hue_col == "Target_Class":
        palette = {"Pseudo-absence": "crimson", "Presence": "dodgerblue"}
    elif palette is None:
        # Dla innych kolumn (np. ciągłych środowiskowych) używamy wbudowanej mapy kolorów
        palette = "viridis"

    plt.figure(figsize=(10, 8))

    # Wykres rozrzutu za pomocą seaborn z dynamicznym parametrem hue
    sns.scatterplot(
        data=df_pca,
        x="PC1",
        y="PC2",
        hue=hue_col,  # Kolorowanie na podstawie podanej zmiennej
        palette=palette,
        alpha=alpha,
        edgecolor="w",
        s=s,
    )

    # Dynamiczne pobieranie procentów wariancji z tablicy wyników
    var_pc1 = explained_variance[0]
    var_pc2 = explained_variance[1]

    # Estetyka i opisy
    plt.title(title, fontsize=14, fontweight="bold", pad=15)
    plt.xlabel(f"Główna składowa 1 (PC1) - {var_pc1:.2f}% wariancji", fontsize=12)
    plt.ylabel(f"Główna składowa 2 (PC2) - {var_pc2:.2f}% wariancji", fontsize=12)

    plt.grid(True, linestyle=":", alpha=0.5)
    plt.legend(title=hue_col, loc="best")
    plt.tight_layout()
    plt.show()



def df_to_map(df, tytul, col_name_to_color):
    """Generuje mapę z punktami pokolorowanymi według wybranej kolumny (klastry lub gatunki).

    df - DataFrame z danymi (musi zawierać 'decimalLongitude' i 'decimalLatitude')
    tytul - string, tytuł wykresu 
    col_name_to_color - string, nazwa kolumny do kolorowania (np. 'Cluster' lub 'species')
    """
    # 1. Inicjalizacja wykresu z odpowiednią projekcją geograficzną
    fig, ax = plt.subplots(
        figsize=(14, 10), subplot_kw={"projection": ccrs.PlateCarree()}
    )

    # Ustawienie szczegółowości mapy
    scale = "50m"

    # Dodajemy naturalne tło lądów, oceanów, granic i wód
    ax.add_feature(cfeature.LAND.with_scale(scale), facecolor="#f5f5f3")
    ax.add_feature(cfeature.OCEAN.with_scale(scale), facecolor="#e0f2fe")
    ax.add_feature(
        cfeature.BORDERS.with_scale(scale),
        edgecolor="#b5b5b5",
        linewidth=0.8,
        linestyle=":",
    )
    ax.add_feature(
        cfeature.LAKES.with_scale(scale),
        facecolor="#e0f2fe",
        edgecolor="#b5b5b5",
        linewidth=0.5,
    )
    ax.add_feature(
        cfeature.RIVERS.with_scale(scale), edgecolor="#cbd5e1", linewidth=0.5
    )

    # =====================================================================
    # AUTOMATYCZNE GENEROWANIE KOLORÓW DLA UNIKALNYCH WARTOŚCI
    # =====================================================================
    unikalne_wartosci = df[col_name_to_color].dropna().unique()

    # Dynamicznie pobieramy paletę kolorów w zależności od liczby grup
    cmap = plt.get_cmap("Set1")

    # 2. Rysujemy punkty dla każdej unikalnej wartości
    for i, wartosc in enumerate(unikalne_wartosci):
        mask = df[col_name_to_color] == wartosc

        # Pobieramy unikalny kolor z palety dla tej grupy
        kolor = cmap(i % cmap.N)

        ax.scatter(
            df.loc[mask, "decimalLongitude"],
            df.loc[mask, "decimalLatitude"],
            c=[kolor],  
            label=str(wartosc),  
            alpha=0.6,
            s=18,
            edgecolors="none",
            transform=ccrs.PlateCarree(),
        )

    # 3. Automatyczne kadrowanie mapy na bazie przekazanego df
    pad = 0.5
    ax.set_extent(
        [
            df["decimalLongitude"].min() - pad,
            df["decimalLongitude"].max() + pad,
            df["decimalLatitude"].min() - pad,
            df["decimalLatitude"].max() + pad,
        ],
        crs=ccrs.PlateCarree(),
    )

    # 4. Dodanie siatki kartograficznej
    gl = ax.gridlines(draw_labels=True, linestyle="--", color="gray", alpha=0.3)
    gl.top_labels = False
    gl.right_labels = False

    # Stylizacja tytułu i legendy
    plt.title(tytul, fontsize=14, fontweight="bold", pad=20)

    plt.legend(
        title=col_name_to_color.capitalize(),
        loc="lower right",
        frameon=True,
        facecolor="white",
        edgecolor="gray",
        framealpha=0.9,
    )

    plt.tight_layout()
    plt.show()


def run_lda_analysis(data_frame, 
                     title="Separacja grup za pomocą analizy LDA", 
                     group_by_col="species",
                     max_components_limit=2, 
                     banned_columns=None):
    """Przeprowadza pełną analizę LDA, wizualizację przestrzeni grup (z procentami na osiach)

    oraz generuje wykresy słupkowe istotności cech środowiskowych.
    """
    # 1. Definiowanie domyślnej listy kolumn do odrzucenia
    if banned_columns is None:
        banned_columns = [
            "gbifID",
            "species",
            "decimalLongitude",
            "decimalLatitude",
            "Cluster",
            "kmeans_cluster",
        ]

    if group_by_col not in banned_columns:
        banned_columns = list(banned_columns) + [group_by_col]

    # Wyodrębnienie cech numerycznych i usunięcie braków danych dla bezpieczeństwa
    df_working = data_frame.dropna().copy()
    features_cols = [col for col in df_working.columns if col not in banned_columns]

    X_numeric = df_working[features_cols].select_dtypes(include=[np.number])
    y_group = df_working[group_by_col]

    if X_numeric.empty:
        print("Błąd: Brak numerycznych cech środowiskowych do analizy.")
        return

    # Standaryzacja danych wejściowych
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_numeric)

    # =====================================================================
    # POPRAWKA: BEZPIECZNE I DYNAMICZNE WYLICZANIE LICZBY SKŁADOWYCH
    # =====================================================================
    num_classes = y_group.nunique()
    num_features = X_numeric.shape[1]
    mathematical_max = min(num_features, num_classes - 1)

    if mathematical_max < 1:
        print(f"Błąd: Zmienna '{group_by_col}' ma tylko {num_classes} klasę/klasy. LDA wymaga minimum 2 klas.")
        return

    # N_comps dopasowuje się do ograniczenia, ale nie przekracza możliwości matematycznych
    n_comps = min(max_components_limit, mathematical_max)

    print("=" * 60)
    print(f"URUCHAMIAM LDA DLA KOLUMNY: '{group_by_col.upper()}'")
    print(f"Liczba unikalnych grup/klas: {num_classes}")
    print(f"Wybrana liczba składowych LDA: {n_comps}")
    print("=" * 60)

    # Uruchomienie modelu LDA
    lda = LinearDiscriminantAnalysis(n_components=n_comps)
    X_lda = lda.fit_transform(X_scaled, y_group)
    
    exp_var = lda.explained_variance_ratio_

    # =====================================================================
    # 3. WYKRES ROZRZUTU (Z PROCENTAMI NA OSIACH)
    # =====================================================================
    plt.figure(figsize=(12, 7))
    for group_name in y_group.unique():
        mask = y_group == group_name
        if n_comps == 2:
            plt.scatter(X_lda[mask, 0], X_lda[mask, 1], label=str(group_name), alpha=0.6, s=25, edgecolors="none")
        else:
            plt.scatter(X_lda[mask, 0], np.zeros_like(X_lda[mask, 0]), label=str(group_name), alpha=0.6, s=25, edgecolors="none")

    # Dynamiczne etykiety osi z dopisanym procentem wyjaśnionej zmienności
    plt.xlabel(f"LDA Składowa 1 ({exp_var[0]*100:.2f}%)", fontsize=11)
    
    if n_comps == 2:
        plt.ylabel(f"LDA Składowa 2 ({exp_var[1]*100:.2f}%)", fontsize=11)
    else:
        plt.gca().get_yaxis().set_visible(False)

    plt.title(title, fontsize=13, fontweight="bold", pad=15)
    plt.legend(title=f"Grupy: {group_by_col}", bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=1.5, fontsize=9)
    plt.tight_layout()
    plt.show()

    # Wyjaśniona zmienność międzyklasowa w konsoli
    print("\n--- WYJAŚNIONA ZMIENNOŚĆ MIĘDZYKLASOWA ---")
    for i, var in enumerate(exp_var):
        print(f"Składowa {i+1} wyjaśnia: {var*100:.2f}% zmienności między grupami")

    # Podpięcie wag cech (Scalings) pod DataFrame
    lda_axes = pd.DataFrame(
        lda.scalings_[:, :n_comps],
        index=X_numeric.columns,
        columns=[f"LDA Składowa {i+1}" for i in range(n_comps)]
    )

    # =====================================================================
    # 4. WYKRESY SŁUPKOWE DLA GŁÓWNYCH SKŁADOWYCH (SCALINGS)
    # =====================================================================
    fig, axes = plt.subplots(1, n_comps, figsize=(7 * n_comps, 8), sharey=False)
    
    if n_comps == 1:
        axes = [axes]
        
    for i, col in enumerate(lda_axes.columns):
        sorted_weights = lda_axes[col].sort_values()
        colors = ['#d9534f' if x < 0 else '#5cb85c' for x in sorted_weights]
        
        axes[i].barh(sorted_weights.index, sorted_weights.values, color=colors, edgecolor='grey', alpha=0.8)
        axes[i].axvline(0, color='black', linestyle='-', linewidth=0.8)
        axes[i].grid(axis='x', linestyle='--', alpha=0.5)
        
        axes[i].set_title(f"Wagi cech dla: {col}\n(wyjaśnia {exp_var[i]*100:.1f}% zmienności)", fontsize=11, fontweight="bold")
        axes[i].set_xlabel("Waga współczynnika (Linear Discriminant Weight)", fontsize=10)
        
    plt.suptitle(f"Wpływ zmiennych środowiskowych na osie LDA ({group_by_col})", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.show()
