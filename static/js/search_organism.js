var scheduled_function = false;
var searchForm = document.getElementById("search-form");
var searchInput = document.getElementById("search-input");
var resultsBox = document.getElementById("results-box");

var csrf = document.getElementsByName("csrfmiddlewaretoken")[0].value;

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
      resultsBox.classList.add("scrollbar");
      if (Array.isArray(data) && searchInput.value.length > 0) {
        resultsBox.classList.add("scrollbar");
        for (var i = 0; i < data.length; i++) {
          var block = document.createElement("button");
          block.setAttribute("id", i);
          block.setAttribute("class", "resultslist");
          block.innerText = data[i]["name"] + " | " + data[i]["class"];
          resultsBox.appendChild(block);
          block.addEventListener("click", function () {
            // alert(this.id)
            let Taxonomy = this.innerText.split(" | ");
            console.log(Taxonomy);
            searchInput.value = Taxonomy[0];
            var setTaxo = document.getElementById("taxo-name");
            setTaxo.value = Taxonomy[0];
            console.log(setTaxo.value);

            resultsBox.innerHTML = "";
            resultsBox.classList.remove("scrollbar");
          });
        }
      } else {
        if (searchInput.value.length > 0) {
          resultsBox.innerHTML = `<b>${data}</b>`;
        } else {
          resultsBox.classList.add("not-visible");
        }
      }
    },
    error: (err) => {
      console.log(err);
    },
  });
}

searchInput.addEventListener("keyup", (e) => {
  resultsBox.innerHTML = "";
  if (resultsBox.classList.contains("not-visible")) {
    resultsBox.classList.remove("not-visible");
  }

  if (scheduled_function) {
    clearTimeout(scheduled_function);
  }

  scheduled_function = setTimeout(function () {
    sendSearchData(e.target.value);
  }, 500);
});
