$(document).ready(function () {
  $(document).on("dblclick", ".editablerow", function () {
    console.log("editable...");
    var value = $(this).text();
    var input =
      "<input type='text' class='input-data' value='" + value + "' />";
    $(this).html(input);
    $(this).removeClass("editablerow");
  });

  $(document).on("blur", ".input-data", function () {
    var value = $(this).val().trim();
    var td = $(this).parent("td");
    $(this).remove();
    td.html(value);
    td.addClass("editablerow");
    var name = td.data("name");
    var type = td.data("type");
    data = { dict_value: name, type: type, value: value };
    console.log(data);
    sendToServer(data, td);
  });
  $(document).on("keypress", ".input-data", function (e) {
    if (e.keyCode === 13) {
      var value = $(this).val().trim();
      var td = $(this).parent("td");
      $(this).remove();
      td.html(value);
      td.addClass("editablerow");
      var name = td.data("name");
      var type = td.data("type");
      data = { dict_value: name, type: type, value: value };
      console.log(data);
      sendToServer(data, td);
    }
  });

  const csrftoken = getCookie("csrftoken");
  const sendToServer = (data, td) => {
    console.log(data);
    $.ajax({
      url: "/dict_update",
      type: "POST",
      headers: { "X-CSRFToken": csrftoken },
      data: data,
    })
      .done((response) => {
        console.log(response);
        if (!response.result) {
          td.append('<p class="text-danger"> You have no permission to change. click here <button class="btn btn-small" onClick="window.location.reload();"> <i class="fa-solid fa-arrows-rotate"></i></button></p>')
        }
      })
      .fail(() => {
        console.log("Error occured");
      });
  };
});
