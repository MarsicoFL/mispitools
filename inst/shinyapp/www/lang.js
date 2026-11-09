var mispiI18n = null;

function mispiI18nTable() {
  if (mispiI18n) return mispiI18n;
  var decode = function(s) {
    var t = document.createElement('textarea');
    t.innerHTML = s;
    return t.value;
  };
  mispiI18n = {};
  i18n_translations.forEach(function(row) {
    var entry = {};
    Object.keys(row).forEach(function(k) {
      if (k !== '_row') entry[k] = decode(row[k]);
    });
    mispiI18n[decode(row._row)] = entry;
  });
  return mispiI18n;
}

function mispiApplyLang(lang) {
  var table = mispiI18nTable();
  $('.i18n').each(function() {
    var row = table[$(this).attr('data-key')];
    if (row && row[lang] !== undefined) $(this).text(row[lang]);
  });
}

$(document).on('click', '.lang-link', function(e) {
  e.preventDefault();
  var lang = $(this).data('lang');
  $('.lang-link').removeClass('active');
  $(this).addClass('active');
  document.documentElement.lang = lang;
  $('input[type=number]').attr('lang', 'en');
  mispiApplyLang(lang);
  Shiny.setInputValue('lang', lang);
});

$(document).on('shiny:connected', function() {
  Shiny.addCustomMessageHandler('i18n-placeholder', function(text) {
    $('.input-group input[type=text][readonly]').attr('placeholder', text);
  });
});
