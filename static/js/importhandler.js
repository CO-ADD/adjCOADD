// custom javascript

$(document).ready(() => {
  console.log("document ready!");
});
// Function to call the function run_task in Views.py
$(".button").on("click", function () {
  $("#preLoader").fadeIn();
  console.log("button clicked");
  const filepath = $("#filepath").text() ? $("#filepath").text() : "none";
  console.log(filepath);
  $.ajax({
    url: "/tasks/",
    data: { type: $(this).data("type"), filepath: filepath.toString() },
    method: "POST",
  })
    .done((res) => {
      console.log(res);
      const html = `
      <tr>
      <td>${res.task_num}</td>
      <td>${res.task_status}</td>
      <td>${res.task_result}</td>
      </tr>`;
      $("#tasks").append(html);
      if (res.status === "error") {
        $("#savedata").prop("disabled", true);
        $("#savedata").addClass("disabled");
      }
      if (res.result == "pass") {
        $("#next_to_confirm").toggleClass("visible");
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
    url: "/tasks/proceed",
    data: {},
    method: "POST",
  })
    .done((res) => {
      $("#preLoader").fadeOut();
      save_data(res);
    })
    .fail((err) => {
      console.log(err);
    });
});

$("#cancelTask").on("click", function () {
  $.ajax({
    url: "/tasks/cancel",
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
  // window.alert("task canceled");
});

save_data = () => {};
