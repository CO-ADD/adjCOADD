var scheduled_function = false;
var searchForm = $("#search-form");
var searchInput = $("#search-input");
var resultsBox = $("#results-box");

var csrf = getCookie('csrftoken')//$("input[name='csrfmiddlewaretoken']").val();

function sendSearchData(inputtext) {
  $.ajax({
    type: "POST",
    url: "/dorganism/search_organism/",
    data: {
      csrfmiddlewaretoken: csrf,
      inputtext: inputtext,
    },
    success: (res) => {
      resultsBox.html("");
      console.log(res);
      const data = res.data;
      resultsBox.addClass("scrollbar");
      if (Array.isArray(data) && searchInput.val().length > 0) {
        resultsBox.addClass("scrollbar");
        for (var i = 0; i < data.length; i++) {
          var block = $("<button>")
            .attr("id", 'taxo_' + i.toString())
            .attr("class", "resultlist")
            .text(data[i]["name"] + " | " + data[i]["class"])
            .appendTo(resultsBox);


          block.on("click", function () {
            console.log("click")
            select_text = $(this).text()
            searchInput.val(select_text);
            var setTaxo = $("#taxo-name");
            let Taxonomy = select_text.split(" | ");
            setTaxo.val(Taxonomy[0]);
            console.log(Taxonomy[0])
            resultsBox.html("");
            resultsBox.removeClass("scrollbar");
          });
        }
      } else {
        if (searchInput.val().length > 0) {
          resultsBox.html(`<b>${data}</b>`);
        } else {
          resultsBox.addClass("not-visible");
        }
      }
    },
    error: (err) => {
      console.log(err);
    },
  });
}

searchInput.on("keyup", (e) => {
  if (searchInput.val().length > 0) {
    if (resultsBox.hasClass("not-visible")) {
      resultsBox.removeClass("not-visible");
    }
    resultsBox.html("start search...<a class='sendtoSearch'>click start</a>");
    // if (scheduled_function) {
    //   clearTimeout(scheduled_function);
    // }
  } else {
    if (resultsBox.hasClass("not-visible") === false) {
      resultsBox.html("")
      resultsBox.addClass("not-visible")
    }
  }

});

$('body').on('click', '.sendtoSearch', () => {
  resultsBox.append('<div class="spinner-border spinner-border-sm" role="status"><span class="visually-hidden">Loading...</span></div>')
  // scheduled_function = setTimeout(function () {
  console.log(searchInput.val())
  sendSearchData(searchInput.val());
  // }, 100);
})