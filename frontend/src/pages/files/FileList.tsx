/**
 * File List Page — TanStack Query migration.
 */
import type React from 'react'
import { useMemo, useState } from 'react'
import { Button, Select, Tag, Modal, Upload, Popconfirm, Tooltip, Space, Typography, Progress } from 'antd'
import {
  UploadOutlined,
  DownloadOutlined,
  FileOutlined,
  InboxOutlined,
  DeleteOutlined,
  FolderOpenOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import type { UploadFile, UploadProps } from 'antd'
import dayjs from 'dayjs'

import {
  useFilesByProject,
  useDeleteFile,
  useBatchUploadFiles,
  useProjectList,
  useSampleList,
} from '@/hooks/queries'
import fileService from '@/services/file.service'
import { PageHeader, DataTable, PageSkeleton, FadeIn, EnhancedEmptyState } from '@/components/common'
import { toast } from '@/utils/notification'
import type { FileItem } from '@/types/file'

const { Option } = Select
const { Title, Text } = Typography

export const FileList: React.FC = () => {
  const [selectedProject, setSelectedProject] = useState<string>('')
  const [selectedSample, setSelectedSample] = useState<string>('')
  const [isUploadModalVisible, setIsUploadModalVisible] = useState(false)
  const [uploadProgress, setUploadProgress] = useState(0)
  const [fileList, setFileList] = useState<UploadFile[]>([])

  // ---- queries -------------------------------------------------------------

  const { data: projectList } = useProjectList({ limit: 200 })
  const projects = (projectList as any)?.items ?? (projectList as any)?.data ?? []

  const { data: filesData, isLoading, isFetching } = useFilesByProject(
    selectedProject || undefined,
  )
  const files: FileItem[] = useMemo(
    () => (filesData as any)?.items ?? (filesData as any) ?? [],
    [filesData],
  )

  const { data: samplesData } = useSampleList(selectedProject || undefined)
  const samples = useMemo(
    () => (samplesData as any)?.items ?? (samplesData as any) ?? [],
    [samplesData],
  )

  // ---- mutations -----------------------------------------------------------

  const deleteMutation = useDeleteFile()
  const uploadMutation = useBatchUploadFiles()

  // ---- handlers ------------------------------------------------------------

  const showUploadModal = () => {
    if (!selectedProject) {
      toast.warning('Please select a project first')
      return
    }
    setIsUploadModalVisible(true)
  }

  const handleCancelUpload = () => {
    setIsUploadModalVisible(false)
    setSelectedSample('')
    setFileList([])
    setUploadProgress(0)
  }

  const uploadProps: UploadProps = {
    name: 'file',
    multiple: true,
    accept: '.fastq,.fq,.fastq.gz,.fq.gz,.bam,.sam,.vcf,.vcf.gz',
    fileList,
    beforeUpload: () => false,
    onChange(info) {
      setFileList(info.fileList)
    },
    onRemove(file) {
      setFileList((prev) => prev.filter((f) => f.uid !== file.uid))
    },
  }

  const handleDelete = async (id: string) => {
    try {
      await deleteMutation.mutateAsync(id)
      toast.success('File deleted successfully')
    } catch {
      toast.error('Failed to delete file')
    }
  }

  const handleDownload = async (file: FileItem) => {
    try {
      await fileService.downloadFile(file.id, file.filename)
    } catch {
      toast.error('Failed to download file')
    }
  }

  const handleUpload = async () => {
    if (!selectedSample) {
      toast.warning('Please select a sample first')
      return
    }
    if (fileList.length === 0) {
      toast.warning('Please select at least one file to upload')
      return
    }
    const filesToUpload: File[] = []
    for (const f of fileList) {
      if (f.originFileObj) filesToUpload.push(f.originFileObj as File)
    }
    if (filesToUpload.length === 0) {
      toast.error('No valid files to upload')
      return
    }
    try {
      await uploadMutation.mutateAsync({
        sampleId: selectedSample,
        files: filesToUpload,
        onProgress: (percent) => setUploadProgress(percent),
      })
      handleCancelUpload()
      toast.success('Files uploaded successfully')
    } catch (error) {
      console.error('Upload error:', error)
    }
  }

  const formatFileSize = (bytes: number) => {
    if (bytes === 0) return '0 B'
    const k = 1024
    const sizes = ['B', 'KB', 'MB', 'GB', 'TB']
    const i = Math.floor(Math.log(bytes) / Math.log(k))
    return `${(bytes / Math.pow(k, i)).toFixed(2)} ${sizes[i]}`
  }

  const columns: ColumnsType<FileItem> = [
    {
      title: 'Filename',
      dataIndex: 'filename',
      key: 'filename',
      render: (name) => (
        <Space>
          <FileOutlined style={{ color: 'var(--color-primary)' }} />
          <span style={{ fontWeight: 500 }}>{name}</span>
        </Space>
      ),
    },
    {
      title: 'Type',
      dataIndex: 'file_type',
      key: 'file_type',
      width: 120,
      render: (type) => {
        const colorMap: Record<string, string> = {
          fastq: 'blue',
          bam: 'green',
          vcf: 'purple',
          sam: 'orange',
        }
        const color = colorMap[type?.toLowerCase()] || 'default'
        return type ? <Tag color={color}>{type.toUpperCase()}</Tag> : '-'
      },
    },
    {
      title: 'Size',
      dataIndex: 'file_size',
      key: 'file_size',
      width: 120,
      render: (size) => formatFileSize(size || 0),
    },
    {
      title: 'Sample',
      dataIndex: 'sample_id',
      key: 'sample_id',
      width: 150,
      render: (sampleId) => {
        const sample = samples.find((s: any) => s.id === sampleId)
        return sample ? sample.sample_id : '-'
      },
    },
    {
      title: 'Uploaded',
      dataIndex: 'uploaded_at',
      key: 'uploaded_at',
      width: 150,
      render: (date) => (date ? dayjs(date).format('YYYY-MM-DD HH:mm') : '-'),
    },
    {
      title: 'Actions',
      key: 'actions',
      width: 150,
      fixed: 'right',
      render: (_, record) => (
        <Space>
          <Tooltip title="Download">
            <Button
              type="text"
              size="small"
              icon={<DownloadOutlined />}
              onClick={() => handleDownload(record)}
            />
          </Tooltip>
          <Popconfirm
            title="Delete File"
            description="Are you sure you want to delete this file?"
            onConfirm={() => handleDelete(record.id)}
            okText="Yes, Delete"
            cancelText="Cancel"
            okButtonProps={{ danger: true, loading: deleteMutation.isPending }}
          >
            <Tooltip title="Delete">
              <Button type="text" size="small" danger icon={<DeleteOutlined />} />
            </Tooltip>
          </Popconfirm>
        </Space>
      ),
    },
  ]

  if (isLoading && files.length === 0 && selectedProject) {
    return <PageSkeleton hasHeader rows={8} />
  }

  return (
    <div>
      <FadeIn direction="up" delay={0} duration={300}>
        <Space direction="vertical" size="middle" style={{ width: '100%', marginBottom: 24 }}>
          <Space align="center">
            <FolderOpenOutlined style={{ fontSize: 28, color: 'var(--color-primary)' }} />
            <div>
              <Title level={3} style={{ margin: 0 }}>
                File Management
              </Title>
              <Text type="secondary">Upload and manage sequencing data files</Text>
            </div>
          </Space>

          <PageHeader
            left={
              <Select
                placeholder="Select Project"
                style={{ width: 300 }}
                value={selectedProject || undefined}
                onChange={setSelectedProject}
                loading={projects.length === 0}
              >
                {projects.map((p: any) => (
                  <Option key={p.id} value={p.id}>
                    {p.name}
                  </Option>
                ))}
              </Select>
            }
            right={
              <Button
                type="primary"
                icon={<UploadOutlined />}
                onClick={showUploadModal}
                disabled={!selectedProject}
              >
                Upload Files
              </Button>
            }
          />
        </Space>
      </FadeIn>

      <FadeIn direction="up" delay={100} duration={300}>
        {!selectedProject ? (
          <EnhancedEmptyState
            type="noData"
            title="Select a project"
            description="Please select a project from the dropdown above to view and manage files"
            size="default"
          />
        ) : files.length === 0 && !isFetching ? (
          <EnhancedEmptyState
            type="noData"
            title="No files uploaded yet"
            description="Upload FASTQ, BAM, SAM, or VCF files to get started with analysis"
            action={{ text: 'Upload Files', onClick: showUploadModal, icon: <UploadOutlined /> }}
            size="default"
          />
        ) : (
          <DataTable
            columns={columns}
            dataSource={files}
            rowKey="id"
            loading={isFetching && files.length === 0}
            pagination={{
              pageSize: 20,
              showSizeChanger: true,
              showTotal: (total) => `Total ${total} files`,
            }}
            emptyText="No Files"
            emptyDescription="Select a project and upload files to get started"
          />
        )}
      </FadeIn>

      <Modal
        title="Upload Files"
        open={isUploadModalVisible}
        onCancel={handleCancelUpload}
        footer={null}
        width={700}
      >
        <div style={{ marginBottom: 16 }}>
          <label style={{ display: 'block', marginBottom: 8, fontWeight: 500 }}>Select Sample *</label>
          <Select
            placeholder="Choose a sample for the files"
            style={{ width: '100%' }}
            value={selectedSample || undefined}
            onChange={setSelectedSample}
          >
            {samples.map((s: any) => (
              <Option key={s.id} value={s.id}>
                {s.sample_id} {s.group_name && `(${s.group_name})`}
              </Option>
            ))}
          </Select>
        </div>

        <Upload.Dragger {...uploadProps} disabled={!selectedSample}>
          <p className="ant-upload-drag-icon">
            <InboxOutlined />
          </p>
          <p className="ant-upload-text">Click or drag files to this area to upload</p>
          <p className="ant-upload-hint">
            Support for single or bulk upload. Accepted formats: FASTQ, BAM, SAM, VCF
            <br />
            Maximum file size: 50GB per file
          </p>
        </Upload.Dragger>

        {uploadMutation.isPending && uploadProgress > 0 && (
          <div style={{ marginTop: 16 }}>
            <Progress percent={uploadProgress} status="active" strokeColor="var(--color-primary)" />
            <Text type="secondary" style={{ fontSize: 12 }}>
              Uploading files... {uploadProgress}%
            </Text>
          </div>
        )}

        <div style={{ marginTop: 16, textAlign: 'right' }}>
          <Space>
            <Button onClick={handleCancelUpload} disabled={uploadMutation.isPending}>
              Cancel
            </Button>
            <Button
              type="primary"
              loading={uploadMutation.isPending}
              disabled={!selectedSample || fileList.length === 0}
              onClick={handleUpload}
            >
              Upload {fileList.length > 0 && `(${fileList.length})`}
            </Button>
          </Space>
        </div>
      </Modal>
    </div>
  )
}

export default FileList
