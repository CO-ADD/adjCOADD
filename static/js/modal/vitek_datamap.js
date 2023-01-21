$(document).ready(function () {
  console.log("loading data vitek");

  $("[data-name=data_process_value]").change(function (e) {
    //selected objects
    var selected_data = [];

    $("input:checkbox[name=type]:checked").each(function () {
      selected_data.push($(this).val().toString());
    });

    // plot value
    var value_str = $(this).val().toString();
    // plot type
    var data_map_str = $("[data-name=data_map_type] option:selected").val()
      .toString;
    // index value
    var index_value = [];
    $("#sortable3 li").each(function () {
      index_value.push($(this).text());
    });
    var index_values_str = index_value.toString();
    // column value
    var column_value = [];
    $("#sortable2 li").each(function () {
      column_value.push($(this).text());
    });
    var column_value_str = column_value.toString();
    // all data
    var data = {
      // data_map: data_map_str,
      selected_data: selected_data,
      values: value_str,
      columns: column_value_str,
      index: index_values_str,
    };
    console.log(data);
    // ajax send data to server, receive data from server
    sendToServer(data, "/vitekcard_list/");
    // var table = result ? result : null;
    // console.log(table);
    // $("#pivotable").html = "";
    // $("#pivotable").html += table;
  });
});

const csrftoken = getCookie("csrftoken");
const sendToServer = (data, url) => {
  console.log("send to server");
  console.log(data);
  var result = "";
  $.ajax({
    url: "/vitekcard_list", //url,
    type: "POST",
    headers: { "X-CSRFToken": csrftoken },
    data: data,
  })
    .done((response) => {
      data = response["table"];
      console.log(data);
      $("#pivotable").html("");
      $("#pivotable").append(data);
    })
    .fail((XMLHttpRequest, textStatus, errorThrown) => {
      console.log(XMLHttpRequest, textStatus, errorThrown);
    });
};
