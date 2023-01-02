// custom javascript
$(document).ready(() => {
  console.log("document ready!");
});
const csrftoken = getCookie("csrftoken");

// Function to call the function run_task in Views.py

$(".button").on("click", function () {
  $("#preLoader").fadeIn();
  console.log("button clicked");
  const filepath = $("#filepath").text() ? $("#filepath").text() : "none";
  const datamodel = $("#datamodel").text() ? $("#datamodel").text() : "none";
  console.log(filepath);
  $.ajax({
    url: "/import/",
    data: {
      type: $(this).data("type"),
      filepath: filepath.toString(),
      datamodel: datamodel.toString(),
    },
    method: "POST",
    headers: { "X-CSRFToken": csrftoken },
  })
    .done((res) => {
      console.log(res);
      // if (jQuery.isEmptyObject(res)) {
      //   console.log("reloading.");
      //   window.location.reload("/import/");
      // }
      const html = `
      <tr>
      <td>${res.task_user}</td>
      <td>${res.task_status}</td>
      <td>${res.task_result}</td>
      </tr>`;
      $("#tasks").append(html);
      if (res.task_status === "Form Errors") {
        $("#save_Proceed").prop("disabled", true);
        $("#save_Proceed").addClass("disabled");
        $("#next_to_confirm").toggleClass("visible");
      }
      if (res.task_status === "Form is Valid") {
        $("#Import_step3").toggleClass("visible");
        $("#progressbar span:first-child").toggleClass("bg-success");
      }
      const taskStatus = res.status;
      $("#preLoader").fadeOut();
    })
    .fail((err) => {
      console.log(err);
    });
});

$("#save_Proceed").on("click", function () {
  $("#preLoader").fadeIn();
  $.ajax({
    url: "/import/",
    method: "POST",
    data: {
      type: $(this).data("type"),
    },
    headers: { "X-CSRFToken": csrftoken },
  })
    .done((res) => {
      $("#preLoader").fadeOut();
      console.log(res);
      // if (res.task_status === "Form is Valid") {
      $("#Import_step4").toggleClass("visible");
      $("#progressbar span:nth-child(2)").toggleClass("bg-success");
      // }
      const html = `<p>${res.status}</p>`;
      $("#mesg_save_Proceed").append(html);
      // save_data(res);
    })
    .fail((err) => {
      console.log(err);
    });
});

$("#stop_Proceed").on("click", function () {
  $.ajax({
    url: "/import/",
    data: {},
    method: "POST",
  })
    .done((res) => {
      window.alert("task canceled");
      location.reload();
    })
    .fail((err) => {
      console.log(err);
    });
});

$("#confirm-save").on("click", function () {
  $.ajax({
    url: "/import/",
    data: {},
    method: "POST",
    data: {
      type: $(this).data("type"),
    },
    headers: { "X-CSRFToken": csrftoken },
  })
    .done((res) => {
      $("#Import_step3").toggleClass("visible");
      $("#Import_step4").toggleClass("visible");
      window.alert(res.status);
      window.location.href = "/import/";
    })
    .fail((err) => {
      console.log(err);
    });
});
