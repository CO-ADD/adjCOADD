var scheduled_function = false;
var searchForm = $("#search-form");
var searchInput = $("#search-input");
var resultsBox = $("#results-box");

var csrf = $("input[name='csrfmiddlewaretoken']").val();

function sendSearchData(inputtext) {
  $.ajax({
    type: "POST",
    url: "/search_organism/",
    data: {
      csrfmiddlewaretoken: csrf,
      inputtext: inputtext,
    },
    success: (res) => {
      console.log(res);
      const data = res.data;
      resultsBox.addClass("scrollbar");
      if (Array.isArray(data) && searchInput.val().length > 0) {
        resultsBox.addClass("scrollbar");
        for (var i = 0; i < data.length; i++) {
          var block = $("<button>")
            .attr("id", i)
            .attr("class", "resultlist")
            .text(data[i]["name"] + " | " + data[i]["class"])
            .appendTo(resultsBox);
          
          block.on("click", function () {
            let Taxonomy = $(this).text().split(" | ");
            console.log(Taxonomy);
            searchInput.val(Taxonomy[0]);
            var setTaxo = $("#taxo-name");
            setTaxo.val(Taxonomy[0]);
            console.log(setTaxo.val());

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

searchInput.on("keyup change", (e) => {
  console.log("keyup event triggered"); 
  resultsBox.html("");
  if (resultsBox.hasClass("not-visible")) {
    resultsBox.removeClass("not-visible");
  }

  if (scheduled_function) {
    clearTimeout(scheduled_function);
  }

  scheduled_function = setTimeout(function () {
    sendSearchData(e.target.value);
  }, 300);
});