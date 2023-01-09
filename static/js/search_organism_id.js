var scheduled_function_id = false;
var searchForm = document.getElementById("id_search-form");
var searchInput = document.getElementById("search-id-input");
var resultsBox = document.getElementById("id_results-box");

var csrf = document.getElementsByName("csrfmiddlewaretoken")[0].value;

function sendSearchData_id(inputtext) {
  $.ajax({
    type: "POST",
    url: "/search_organism_id/",
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
        for (let i = 0; i < data.length; i++) {
          let block = document.createElement("button");
          block.setAttribute("id", i);
          block.setAttribute("class", "resultslist");
          block.innerText = data[i]["name"];
          resultsBox.appendChild(block);
          block.addEventListener("click", function () {
            let organism_id = this.innerText;
            searchInput.value = organism_id;
            let setMicro = document.getElementById("microorg");
            setMicro.value = organism_id;
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
  console.log("searching....");
  resultsBox.innerHTML = "";
  if (resultsBox.classList.contains("not-visible")) {
    resultsBox.classList.remove("not-visible");
  }
  // sendSearchData(e.target.value);

  if (scheduled_function_id) {
    clearTimeout(scheduled_function_id);
  }

  scheduled_function_id = setTimeout(function () {
    sendSearchData_id(e.target.value);
  }, 500);
});
