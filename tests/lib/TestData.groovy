import java.nio.file.*
import java.net.*

class TestData {

  static File stageDir(String launchDir, String relative) {
    def d = new File(launchDir, relative)
    d.mkdirs()
    return d
  }

  static void fetch(String url, File dst, Map opts = [:]) {
    dst.parentFile?.mkdirs()
    def conn = new URL(url).openConnection()
    conn.setRequestProperty("User-Agent", opts.ua ?: "nf-test")
    conn.connectTimeout = (opts.timeout ?: 15000) as int
    conn.readTimeout    = (opts.timeout ?: 15000) as int
    dst.withOutputStream { os ->
      conn.inputStream.withCloseable { it.transferTo(os) }
    }
  }

  static void fetchAll(String base, List<String> files, File destDir) {
    files.each { fname -> fetch("${base}/${fname}", new File(destDir, fname)) }
  }


  static void fetchCached(String base, List<String> files, File destDir, File cacheDir) {
    cacheDir.mkdirs()
    files.each { fn ->
      def cached = new File(cacheDir, fn)
      if (!cached.exists()) fetch("${base}/${fn}", cached)
      Files.copy(cached.toPath(), new File(destDir, fn).toPath(),
                 StandardCopyOption.REPLACE_EXISTING)
    }
  }
}
