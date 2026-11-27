# Українська редакція рукопису AAM

Робочий переклад для спільного наукового обговорення з професором Богдановим.
Джерело — `paper/submissions/aam/manuscript.tex` у коміті
`7597cba1a378e6f436a3bfab304f85702902ee4e` репозиторію CrackFromHoleLEFM.
Відповідність Git blob, SHA-256, історію PR і відмінності пізнішої гілки
зафіксовано в [аудиті джерела](audit/SOURCE_REVIEW.md).

[Попередній PDF PR #89 (НЕ відповідає новому TeX)](AAM_UA_Bogdanov_review.pdf) ·
[Початкові файли та ідентифікатори](SOURCE_MANIFEST.json) ·
[Наукові питання](audit/SCIENTIFIC_QUESTIONS.md)

## Відтворення PDF

Потрібні XeLaTeX, BibTeX, STIX Two Text (Regular, Bold, Italic, Bold Italic), STIX Two Math, пакет unicode-math та українські шаблони переносу слів. Файли STIXTwoText-*.otf і STIXTwoMath-Regular.otf мають бути доступні XeLaTeX (наприклад, пакет stix2-otf у TeX Live). Noto Serif залишається потрібним лише для повторної локалізації рисунків. З кореня репозиторію у PowerShell:

```powershell
.\paper\ukrainian_review\build_ua.ps1
```

Сценарій виконує XeLaTeX → BibTeX → XeLaTeX → XeLaTeX і перевіряє
журнал перед копіюванням PDF. Еквівалентна ручна збірка з каталогу
`paper/ukrainian_review`:

```text
xelatex -halt-on-error -recorder -jobname=AAM_UA_Bogdanov_review main_ua.tex
bibtex AAM_UA_Bogdanov_review
xelatex -halt-on-error -recorder -jobname=AAM_UA_Bogdanov_review main_ua.tex
xelatex -halt-on-error -recorder -jobname=AAM_UA_Bogdanov_review main_ua.tex
```

Файли `main_ua.tex`, `references.bib` і `figures/` достатні для
збірки за наявності вказаних залежностей. Сценарії аудиту використовують
Git-об’єкти оригіналу, тому для повного аудиту потрібна копія
репозиторію з історією, яка містить базовий коміт.

## Рисунки та перевірки

Збірка включає 31 підготовлений векторний PDF. Відтворення українських
написів потребує Python з PyMuPDF та pypdf, оригінальних PDF з `paper/figures`
та Noto Serif Regular TTF. Приклад з кореня репозиторію:

```text
python paper/ukrainian_review/scripts/localize_figures.py --font NotoSerif-Regular.ttf
python paper/ukrainian_review/scripts/verify_preservation.py --write-audit
python paper/ukrainian_review/scripts/verify_review_pdf.py paper/ukrainian_review/.build --write-audit
```

Перший сценарій працює із затвердженими векторними експортами та
замінює лише текстові позначення. Сценарій зберігає затверджені
геометрію, осі, маркери, криві, дані сіток і десяткову точність
чисел. Для відновлення з тими самими розмірами тексту використовуйте
версії залежностей та ідентичний шрифт із
[паспорта рисунків](audit/FIGURE_PROVENANCE.json).
Скінченно-елементний розрахунок та MATLAB для збірки не потрібні.

Етап перекладу PR #89 завершено на `REVIEW_STOP_AAM_UA_TRANSLATION`. За вказівкою дослідника авторами робочого рукопису зазначено Вячеслава Богданова та Данила Селіванова у цьому порядку. Поточна редакторська гілка містить відредагований TeX, але **повну повторну компіляцію ще не виконано**. PDF та BUILD_VERIFICATION.json успадковані з PR #89 і є історичними доказами, а не результатами перевірки нового джерела. Перед переданням Claude або Богданову необхідно зібрати і перевірити новий PDF. Див. [редакторський аудит](audit/PRE_CLAUDE_REVIEW.md).
